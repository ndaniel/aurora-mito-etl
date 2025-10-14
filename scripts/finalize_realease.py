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
import requests
import math
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
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

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import pandas as pd
import numpy as np
import math

def add_tanimoto_scores(
    df: pd.DataFrame,
    fp_radius: int = 2,      # 2 → ECFP4
    nbits: int = 2048,
    topk: int = 3
) -> pd.DataFrame:
    """
    Append RDKit-based similarity scores + confidence label to a mito.txt-like DataFrame.

    Input df must have columns: 'compound', 'known_status', 'SMILES'
    Adds columns:
      - MaxSim_all (float or None)
      - TopKMean_all (float or None)
      - BestRef_name (str or None)
      - confidence_similarity (str in {'high','medium','low','very-low'} or None)
    """

    # ----------------- helpers -----------------
    def mol_from_smiles(s):
        if not isinstance(s, str) or not s.strip():
            return None
        m = Chem.MolFromSmiles(s)
        if m is None:
            return None
        Chem.SanitizeMol(m)
        return m

    def ecfp(m):
        return AllChem.GetMorganFingerprintAsBitVect(m, radius=fp_radius, nBits=nbits)

    def inchikey(m):
        try:
            return Chem.MolToInchiKey(m)
        except Exception:
            return None

    def topk_mean(vals, k=3):
        if not vals:
            return np.nan
        arr = sorted(vals, reverse=True)
        k = min(k, len(arr))
        return float(np.mean(arr[:k]))

    def confidence_from_sims(maxsim, topkmean):
        # Prefer TopK mean; fall back to MaxSim if TopK is NaN
        v = topkmean if (topkmean == topkmean) else maxsim  # NaN check
        if v != v or v is None:          # both NaN/None
            return None
        if v >= 0.70:
            return "high"
        if v >= 0.50:
            return "medium"
        if v >= 0.30:
            return "low"
        return "very-low"

    # ----------------- prepare refs -----------------
    refs = df[(df["known_status"].astype(str).str.lower() == "known") & df["SMILES"].notna()].copy()
    refs["mol"] = refs["SMILES"].apply(mol_from_smiles)
    refs = refs.dropna(subset=["mol"]).copy()
    refs["inchikey"] = refs["mol"].apply(inchikey)
    refs = refs.drop_duplicates(subset=["inchikey"]).reset_index(drop=True)
    refs["fp"] = refs["mol"].apply(ecfp)

    # If no valid refs, just append empty columns and return
    if refs.empty:
        out = df.copy()
        out["MaxSim_all"] = None
        out["TopKMean_all"] = None
        out["BestRef_name"] = None
        out["confidence_similarity"] = None
        return out

    ref_fps   = refs["fp"].tolist()
    ref_names = refs["compound"].tolist()

    # ----------------- compute per-row -----------------
    rows = []
    for _, row in df.iterrows():
        smi = row.get("SMILES", None)
        max_sim = np.nan
        mean_tk = np.nan
        best_name = None

        if isinstance(smi, str) and smi.strip():
            mol = mol_from_smiles(smi)
            if mol is not None:
                fp = ecfp(mol)
                sims = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
                if sims:
                    max_idx = int(np.argmax(sims))
                    max_sim = float(sims[max_idx])
                    mean_tk = topk_mean(sims, topk)
                    best_name = ref_names[max_idx]

        rows.append({
            "MaxSim_all": round(max_sim, 3) if (max_sim == max_sim) else None,
            "TopKMean_all": round(mean_tk, 3) if (mean_tk == mean_tk) else None,
            "BestRef_name": best_name,
            "confidence_similarity": confidence_from_sims(max_sim, mean_tk),
        })

    return pd.concat([df.reset_index(drop=True), pd.DataFrame(rows)], axis=1)


# --------------------------------------------------------------------------
# Define output file path
# --------------------------------------------------------------------------
OUTPUT = release_dir / "all_mito_complex_I_inhibitors.txt"
OUTPUT_EXCEL = release_dir / "all_mito_complex_I_inhibitors.xlsx"
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
    norm = norm.replace("\u2013", "-").replace("\u2014", "-").replace("\u2212", "-") # minus sign
    norm = re.sub(r"\s+", " ", norm)               
    # strip trailing 's' if plural and not an obvious chemical suffix (e.g., 'us', 'is')
    if len(norm) > 4 and norm.endswith("s") and not norm.endswith(("us", "is","os","gas")):
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

ref = set(unique_refs.values())
ref.difference_update(["Roterone","Piericidin","Bongkrekic"])
ref.add("Piericidin A")
ref.add("Bongkrekic acid")
ref = sorted(ref)

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
    blacklist=set(["zinc","complex i","complex 1","complex i blockers","complex i blocker","complex i inhibitor","complex i inhibitors","iron","compound","components of cigarette smoke","hepatitis c virus","roterone","ozone","calcium","no small-molecule compound named","SiO2 nanoparticles","crude oil","derivatives","dispersed oil","extensively oxidized low-density lipoprotein","fatty acids","hydrogen gas","inhibitor derived from nadh","inorganic arsenic","lithium","methane","nickel ion","nitrate","vehicle of sandimmun","camel milk exosomes","acidic buffer","mir-27a-3p","fish oil","cadmium","arsenic trioxide","chromium","hexavalent chromium"])
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

# add smiles
###############################################################################


targets = stats.loc[stats["confidence"]!="low","compound"].unique().tolist()

def _first_smiles(props: dict):
    """Pick the best available SMILES key from a PubChem properties dict."""
    for k in ("IsomericSMILES", "CanonicalSMILES", "SMILES", "ConnectivitySMILES"):
        if k in props and props[k]:
            return props[k], k
    return None, None

def fetch_smiles(name: str, normalize : bool = True):
    """Try PubChem, then fallback to ChEMBL (URL-safe, with spaces encoded)."""
    if normalize:
        base = normalize_compound_name(name)
    else:
        base = name.strip().replace("\u2013", "-").replace("\u2014", "-").replace("\u2212", "-") # minus sign
    encoded = requests.utils.quote(base, safe="")  # ensure spaces → %20

    # --- PubChem ---
    try:
        url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/property/IsomericSMILES,CanonicalSMILES,SMILES,ConnectivitySMILES,InChIKey/JSON")
        r = requests.get(url, timeout=20)
        if r.ok:
            props = r.json()["PropertyTable"]["Properties"][0]
            smi, which = _first_smiles(props)
            if smi:
                return smi, f"PubChem:{which}", base
    except Exception:
        pass

    # --- ChEMBL fallback ---
    try:
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search?q={encoded}&format=json"
        r = requests.get(url, timeout=20)
        if r.status_code == 200:
            data = r.json()
            if data.get("molecules"):
                smi = data["molecules"][0].get("molecule_structures", {}).get("canonical_smiles")
                if smi:
                    return smi, "ChEMBL", base
    except Exception:
        pass

    return None, None, base


smiles = {}
i = 0
n = len(targets)
for t in targets:
    i = i + 1
    smi, source, query = fetch_smiles(t)
    if not smi:
        # try one more time
        smi, source, query = fetch_smiles(t,normalize=False)          
    print(f"{i}/{n} compounds: {t} -> {smi} - {source} - {query}")  
    smiles[t] = smi
    time.sleep(0.2)

stats["SMILES"] = stats["compound"].map(smiles).fillna("")

#
stats = add_tanimoto_scores(stats)

# Keep SMILES as the final column in the exported table for readability
if "SMILES" in stats.columns:
    ordered_cols = [c for c in stats.columns if c != "SMILES"] + ["SMILES"]
    stats = stats[ordered_cols]

# Write stats to OUTPUT
stats.to_csv(OUTPUT, sep="\t", index=False)
stats.to_excel(OUTPUT_EXCEL, index=False)

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

append_release_info(
    OUTPUT_EXCEL,
    step=step_name,
    sources=common_sources,
    parameters={
        **common_params,
        "known_ref_boost": 100,
        "confidence_bins": "low: ≤1; low-medium: 2; medium: 3–5; high: ≥6",
        "pubmed_ids_concat": "unique PMIDs, sorted, ';' separated",
    },
    notes=(
        "Excel mirror of the aggregated compound summary with SMILES and RDKit "
        "similarity metrics."
    ),
)

print('[INFO] Setup complete.')
