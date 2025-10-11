#!/usr/bin/env python3
"""
Batched + resumable GPT-4.1-mini classification of PubMed abstracts.

Input  (TSV): data/staging/pubmed_filtered.txt
  columns: pmid<TAB>year<TAB>title_and_abstract

Output (TSV): data/staging/pubmed_gpt.txt
  columns: PMID<TAB>YES|probablyYES|NO<TAB>compound_names

Behavior:
- Batches abstracts into single API calls.
- Strict TSV-only output (no explanations).
- Normalizes/forces PMID format and deduplicates.
- Resumable via reading existing output (can be disabled with FORCE_RERUN).
"""

import csv
import os
import re
import time
from pathlib import Path

import openai
from tqdm import tqdm

# ------------------------------------------------------------------------------
# Repo-aware paths (etl/... -> repo root)
# ------------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
INPUT  = ROOT / "data" / "staging" / "pubmed_filtered.txt"
OUTPUT = ROOT / "data" / "staging" / "pubmed_gpt.txt"

# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------
MODEL = "gpt-4.1-mini"     # cheapest 4.x mini model
BATCH_SIZE = 10            # abstracts per API request
SLEEP_BETWEEN = 1.0        # seconds between successful batches
MAX_TOKENS_PER_ITEM = 25   # reply budget per abstract (short TSV lines)
FORCE_RERUN = False        # set True to start from 0 even if output exists
TEXT_MAX_CHARS = 6000      # truncate very long abstracts for token safety

openai.api_key = os.getenv("OPENAI_API_KEY")
assert openai.api_key, "Missing OPENAI_API_KEY"

# ------------------------------------------------------------------------------
# Utilities
# ------------------------------------------------------------------------------

BANNED_SIMPLE = {
    "na","k","cl","mg","ca","fe","cu","zn","mn","au","ag","h2o","o2"
    # add "no" here if you want to forbid nitric oxide as a compound output
}

# Accept 1–9 digits; allow optional "PMID" prefix; anchor to the whole field
PMID_RX_STRICT = re.compile(r'^(?:PMID)?\s*(\d{1,9})$')
# Fallback: find a 1–9 digit chunk anywhere (used only if strict fails)
PMID_RX_LOOSE  = re.compile(r'(\d{1,9})')

def norm_pmid(s: str) -> str:
    """Extract a PMID (1–9 digits). Prefer exact field match; fall back to search."""
    if not s:
        return ""
    m = PMID_RX_STRICT.match(s)
    if m:
        return m.group(1)
    m = PMID_RX_LOOSE.search(s)
    return m.group(1) if m else ""


def clean_text(s: str) -> str:
    """Normalize whitespace and truncate to keep token use sane."""
    if not s:
        return ""
    s = " ".join(s.split())  # collapse whitespace, including newlines
    if len(s) > TEXT_MAX_CHARS:
        s = s[:TEXT_MAX_CHARS]
    return s

def load_done_pmids(path: Path) -> set[str]:
    """Load already processed PMIDs from an existing output file (normalized)."""
    done: set[str] = set()
    if not path.exists() or FORCE_RERUN:
        return done
    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            first_field = line.split("\t", 1)[0]
            pmid = norm_pmid(first_field)
            if pmid:
                done.add(pmid)
    return done

# ------------------------------------------------------------------------------
# Model I/O
# ------------------------------------------------------------------------------
def ask_batch(batch: list[tuple[str, str]]) -> list[str]:
    """
    Ask the model about multiple abstracts at once.
    Returns exactly len(batch) lines of TSV:
      PMID<digits>\tYES|probablyYES|NO\tcompound_names|NA
    """
    joined = "\n\n".join(
        f"[{i+1}] PMID {pmid}\n{text}"
        for i, (pmid, text) in enumerate(batch)
    )

    system_msg = (
        "Extract claims from biomedical abstracts. Decide ONLY whether a SMALL-MOLECULE "
        "compound inhibits mitochondrial Complex I (respiratory Complex I, CI, NADH "
        "dehydrogenase, NADH:ubiquinone oxidoreductase). "
        "Ignore peptides/proteins/antibodies, gene silencing, mutations, polymers, and other complexes (II–V). "
        "Treat clear decreases of Complex I activity or mitochondrial NADH oxidation as probablyYES. "
        "Never output elements/ions/simple salts (Na, K, Cl, Mg, Ca, Fe, Cu, Zn, Mn, Au, Ag, etc.). "
        "Output must be TSV lines only—no explanations or extra text."
    )

    user_prompt = (
        "For EACH numbered abstract below, output EXACTLY one line, in the SAME ORDER. "
        "Each line MUST start with the literal string 'PMID' immediately followed by the digits of that PMID, "
        "then a TAB, then ONE of {YES, probablyYES, NO}, then a TAB, then compound names "
        "(up to 3, separated by '; '), or 'NA' if unclear or NO.\n"
        "Format:\n"
        "PMID<digits>\tYES|probablyYES|NO\t<compound_names|NA>\n"
        "Examples:\n"
        "PMID12345678\tYES\trotenone\n"
        "PMID23456789\tprobablyYES\tmetformin\n"
        "PMID34567890\tYES\tpiericidin A; fenpyroximate\n"
        "PMID45678901\tNO\tNA\n\n"
        f"{joined}"
    )

    resp = openai.chat.completions.create(
        model=MODEL,
        messages=[
            {"role": "system", "content": system_msg},
            {"role": "user",  "content": user_prompt},
        ],
        temperature=0,
        max_tokens=MAX_TOKENS_PER_ITEM * len(batch),
    )
    return resp.choices[0].message.content.strip().splitlines()

def parse_answer(line: str, expected_pmid: str) -> tuple[str, str, str] | None:
    """
    Parse one TSV line into (pmid, flag, compounds).
    Aligns/forces PMID to expected PMID and sanitizes compounds.
    """
    fields = line.strip().split("\t")
    if len(fields) < 3:
        return None

    out_pmid_raw, flag_raw, compounds_raw = (f.strip() for f in fields[:3])

    # Normalize PMID and force to expected if the model drifts.
    out_pmid = norm_pmid(out_pmid_raw) or expected_pmid
    if out_pmid != expected_pmid:
        out_pmid = expected_pmid

    flag_l = flag_raw.strip().lower()
    if flag_l not in {"yes", "probablyyes", "no"}:
        return None

    # sanitize compounds (max 3; strip junk; drop simple ions)
    names: list[str] = []
    for name in (c.strip() for c in compounds_raw.split(";")):
        if not name:
            continue
        nlow = name.lower()
        if nlow in BANNED_SIMPLE:
            continue
        if len(name) > 80:  # guard against explanations creeping in
            continue
        name = name.strip(" ,.;:()[]{}")
        if not name:
            continue
        names.append(name)
        if len(names) == 3:
            break

    if flag_l == "no":
        return out_pmid, "NO", "NA"
    if flag_l == "yes":
        return out_pmid, "YES", ("; ".join(names) if names else "NA")
    return out_pmid, "probablyYES", ("; ".join(names) if names else "NA")

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
def main() -> None:
    print(f"[INFO] Starting GPT-4.1-mini classification (batch={BATCH_SIZE})")
    print(f"[INFO] Input : {INPUT}")
    print(f"[INFO] Output: {OUTPUT}")

    done_pmids = load_done_pmids(OUTPUT)
    if done_pmids and not FORCE_RERUN:
        print(f"[INFO] Resuming — {len(done_pmids):,} PMIDs already done")
    else:
        if FORCE_RERUN and OUTPUT.exists():
            print("[INFO] FORCE_RERUN=True — existing output will be appended to, but no PMIDs will be skipped based on it")

    with INPUT.open(newline="", encoding="utf-8") as fin, \
         OUTPUT.open("a", newline="", encoding="utf-8") as fout:

        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")

        batch: list[tuple[str, str]] = []
        total_sent = 0
        backoff = 10

        for row in tqdm(reader, desc="Scanning input"):
            if len(row) < 3:
                continue
            pmid = norm_pmid(row[0])
            if not pmid:
                continue
            if pmid in done_pmids:
                continue

            text = clean_text(row[2])
            if not text:
                continue

            batch.append((pmid, text))

            if len(batch) >= BATCH_SIZE:
                try:
                    answers = ask_batch(batch)

                    # Align answers to inputs; dedupe within this run and across runs
                    seen_this_run: set[str] = set()
                    for (pmid_expected, _text), line in zip(batch, answers):
                        parsed = parse_answer(line, expected_pmid=pmid_expected)
                        if not parsed:
                            continue
                        pmid_norm = parsed[0]
                        if pmid_norm in seen_this_run or pmid_norm in done_pmids:
                            continue
                        writer.writerow(parsed)
                        seen_this_run.add(pmid_norm)
                        done_pmids.add(pmid_norm)

                    fout.flush()
                    total_sent += len(batch)
                    batch = []
                    time.sleep(SLEEP_BETWEEN)
                    backoff = 10
                except Exception as e:
                    print(f"[WARN] API error: {e}; sleeping {backoff}s…")
                    time.sleep(backoff)
                    backoff = min(backoff * 2, 120)

        # last partial batch
        if batch:
            answers = ask_batch(batch)
            seen_this_run: set[str] = set()
            for (pmid_expected, _text), line in zip(batch, answers):
                parsed = parse_answer(line, expected_pmid=pmid_expected)
                if not parsed:
                    continue
                pmid_norm = parsed[0]
                if pmid_norm in seen_this_run or pmid_norm in done_pmids:
                    continue
                writer.writerow(parsed)
                seen_this_run.add(pmid_norm)
                done_pmids.add(pmid_norm)
            fout.flush()
            total_sent += len(batch)

    print(f"[DONE] Processed ~{total_sent:,} new abstracts; total {len(done_pmids):,} lines → {OUTPUT}")

if __name__ == "__main__":
    main()

