#!/usr/bin/env python3
"""
Finalize processed output for mitochondrial complex I inhibitors.

- Creates a versioned folder under data/processed/<YYYY-MM-DD>/
- Deletes old files under data/processed/
- Defines canonical input/output paths for later processing steps
- Appends provenance entries to data/processed/<date>/release_info.txt
"""

import os
import shutil
import time
import re
import hashlib
import datetime as _dt
import pandas as pd
import numpy as np
from pathlib import Path

# --------------------------------------------------------------------------
# Repo-aware paths
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
STAGING_GPT = ROOT / "data" / "staging" / "pubmed_gpt.txt"
REF_INHIBITORS = ROOT / "data" / "reference" / "mitochondrial_complex_i_inhibitors.txt"
PROCESSED_DIR = ROOT / "data" / "processed"

# --------------------------------------------------------------------------
# Helpers: checksum + provenance appender
# --------------------------------------------------------------------------
def sha256sum(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()

def append_release_info(artifact_path: Path, step: str, sources: list[str], parameters: dict | None = None, notes: str | None = None) -> None:
    """
    Append a YAML-ish block describing an artifact into a release_info.txt
    located in the artifact's parent directory.
    """
    info_path = artifact_path.parent / "release_info.txt"
    checksum = sha256sum(artifact_path)
    today = _dt.date.today().isoformat()

    # Render sources and parameters lists
    src_lines = "\n".join([f"  - {s}" for s in sources]) if sources else "  - n/a"
    if parameters:
        param_lines = "\n".join([f"  {k}: {v}" for k, v in parameters.items()])
    else:
        param_lines = "  {}"

    block = (
        "# ETL artifact provenance\n"
        f"file: {artifact_path.name}\n"
        f"sha256: {checksum}\n"
        f"generated_on: {today}\n"
        f"step: {step}\n"
        "sources:\n"
        f"{src_lines}\n"
        "parameters:\n"
        f"{param_lines}\n"
    )
    if notes:
        block += "notes: |\n" + "\n".join([f"  {line}" for line in notes.splitlines()]) + "\n"
    block += "\n"

    # Append (create file if missing)
    with open(info_path, "a", encoding="utf-8") as f:
        f.write(block)

# --------------------------------------------------------------------------
# Create dated processed folder and clean old outputs
# --------------------------------------------------------------------------
today = time.strftime("%Y-%m-%d")
release_dir = PROCESSED_DIR / today

print(f"[INFO] Removing old processed files from: {PROCESSED_DIR}")
for item in PROCESSED_DIR.glob("*"):
    if item.is_file():
        item.unlink()
    elif item.is_dir():
        shutil.rmtree(item)

print(f"[INFO] Creating release folder: {release_dir}")
release_dir.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# Define output file path
# --------------------------------------------------------------------------
OUTPUT = release_dir / "all_mito_complex_I_inhibitors.txt"
OUTPUT_NEW = release_dir / "new_mito_complex_I_inhibitors.txt"
OUTPUT_NEW_EXCEL = release_dir / "new_mito_complex_I_inhibitors.xlsx"

print(f"[INFO] Output will be written to: {OUTPUT}")

print("[INFO] Reading known complex I inhibitors.")
ref = []
with open(REF_INHIBITORS, "r", encoding="utf-8") as f:
    ref = sorted(set([e.strip() for e in f if e.strip()]))

# --------------------------------------------------------------------------
# Deduplicate known reference compounds (normalize: remove spaces/dashes, lowercase)
# --------------------------------------------------------------------------
def normalize_compound_name(name: str) -> str:
    # remove spaces and dashes, lowercase
    norm = re.sub(r"[-\s]", "", name.strip().lower())
    # strip trailing 's' if plural and not an obvious chemical suffix (e.g., 'us', 'is')
    if len(norm) > 4 and norm.endswith("s") and not norm.endswith(("us", "is")):
        norm = norm[:-1]
    return norm

unique_refs = {}
for c in ref:
    if not c.strip():
        continue
    norm = normalize_compound_name(c)
    # keep the first version seen
    if norm not in unique_refs:
        unique_refs[norm] = c.strip()

ref = sorted(unique_refs.values())

def openparanthese(s: str) -> str:
    #
    ix1 = s.find("(")
    ix2 = s.find(")")
    r = s
    if ix1!=-1 and ix2==-1:
        r = s.partition("(")[0]
    return r

print("[INFO] Reading newly found complex I inhibitors from PubMed.")
inh = []
with open(STAGING_GPT, "r", encoding="utf-8") as f:
    inh = [e.rstrip("\r\n").split("\t") for e in f if e.rstrip("\r\n")]
    inh = [e for e in inh if e[1] and e[1].lower() != "no" and e[2] and e[2].lower() != "na"]
    # split the ; on several rows
    r = []
    for e in inh:
        if e[2].find(";")==-1:
            r.append(e)
        else:
            t = [[e[0],e[1],x.strip()] for x in e[2].split(";") if x.strip()]
            r.extend(t)
    inh = r[:]
    inh = [(e[0],e[1],openparanthese(e[2])) for e in inh if e[2].strip()]
    inh = [(e[0],e[1],e[2].replace("analogs","").replace("analogue","").replace("analog","").replace("diphenyleneiodonium","diphenylene iodonium").replace("acetogenins","acetogenin")) for e in inh if e[2].strip()]
    inh = [(e[0],e[1],e[2].strip()) for e in inh if e[2].strip()]
    inh = [e for e in inh if e[1] and e[1].lower() != "no" and e[2] and e[2].lower() != "na"]
    blacklist=set(["zinc","complex i","complex 1","complex i blockers","complex i blocker","complex i inhibitor","complex i inhibitors","iron","compound","components of cigarette smoke","hepatitis c virus"])
    inh = [e for e in inh if len(e[2]) > 2 and e[2].lower() not in blacklist]
    inh = [e for e in inh if e[2].lower()!="no" and e[2].lower().find("nitric oxide")==-1 and e[2].lower().find("mitochondr")==-1 and e[2].lower().find("silencing")==-1 and not e[2].lower().startswith("compound")]

# some statistics
df = pd.DataFrame(inh, columns=["pmid", "confidence", "compound"])

# Basic cleaning
df["pmid"] = df["pmid"].astype(str).str.strip()
df["confidence"] = df["confidence"].str.strip()
df["compound"] = df["compound"].str.strip()

# Write stats to OUTPUT_NEW
df.to_csv(OUTPUT_NEW, sep="\t", index=False)

# Write stats to OUTPUT_NEW_EXCEL (clickable hyperlinks)
dfe = df.copy()
dfe["URL"] = '=HYPERLINK("https://pubmed.ncbi.nlm.nih.gov/' + dfe["pmid"].astype(str) + '/","' + dfe["pmid"].astype(str) + '")'
dfe.to_excel(OUTPUT_NEW_EXCEL, index=False)

# Build stats: compound, pubmed_references, known_status
# - Count UNIQUE PMIDs per normalized compound
df["_compound_norm"] = df["compound"].str.lower()

stats = (
    df.groupby("_compound_norm", as_index=False)
      .agg(
          pubmed_references=("pmid", "nunique"),
          compound=("compound", "first"),
          pubmed_ids=("pmid", lambda x: ";".join(map(str, sorted(set(x)))))
      )[["compound", "pubmed_references", "pubmed_ids"]]
      .sort_values(["pubmed_references", "compound"], ascending=[False, True])
      .reset_index(drop=True)
)

# Per your instruction: set all to "new"
stats["known_status"] = "new"

# --- Append all known reference compounds ---
known_df = pd.DataFrame({
    "compound": ref,
    "pubmed_references": [100] * len(ref),
    "pubmed_ids": "",
    "known_status": ["known"] * len(ref),
})

# Concatenate both
stats = pd.concat([stats, known_df], ignore_index=True)

stats["confidence"] = pd.cut(
    stats["pubmed_references"],
    bins=[-np.inf, 1, 2, 5, np.inf],
    labels=["low", "low-medium", "medium", "high"]
)

stats = stats[["compound", "pubmed_references", "known_status", "confidence", "pubmed_ids"]]


# Final sort
stats = stats.sort_values(["pubmed_references", "compound"], ascending=[False, True]).reset_index(drop=True)

# Write stats to OUTPUT
stats.to_csv(OUTPUT, sep="\t", index=False)

print("[INFO] Writing release_info entries.")
step_name = Path(__file__).name
common_sources = [
    f"LLM classifications: {STAGING_GPT.relative_to(ROOT)}",
    f"Known inhibitors list: {REF_INHIBITORS.relative_to(ROOT)}",
]
common_params = {
    "script": step_name,
    "date": today,
}

# Append provenance for each artifact
append_release_info(
    OUTPUT_NEW,
    step=step_name,
    sources=common_sources,
    parameters={**common_params},
    notes=(
        "Newly found inhibitors with per-PMID records.\n"
        "Columns: pmid, confidence (LLM label), compound.\n"
        "Excel version also includes a clickable PubMed URL."
    ),
)

append_release_info(
    OUTPUT_NEW_EXCEL,
    step=step_name,
    sources=common_sources,
    parameters={**common_params},
    notes=(
        "Excel mirror of new inhibitors with clickable PubMed hyperlinks "
        'in column \"URL\".'
    ),
)

append_release_info(
    OUTPUT,
    step=step_name,
    sources=common_sources,
    parameters={
        **common_params,
        "known_ref_boost": 100,
        "confidence_bins": "low: ≤1; low-medium: 2; medium: 3–5; high: ≥6",
        "pubmed_ids_concat": "unique PMIDs, sorted, ';' separated",
    },
    notes=(
        "Aggregated stats per compound including:\n"
        "- pubmed_references (unique PMID count)\n"
        "- known_status (reference list appended at 100)\n"
        "- confidence: categorical score based on pubmed_references count "
        "(low, low-medium, medium, high)\n"
        "- pubmed_ids (unique, sorted, ';' separated)"
    ),
)

print('[INFO] Setup complete.')

