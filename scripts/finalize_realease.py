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
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, inchi
from rdkit.Chem import rdFingerprintGenerator as rfg
from rdkit.Chem.rdmolops import RemoveHs



# --------------------------------------------------------------------------
# Repo-aware paths
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
STAGING_GPT = ROOT / "data" / "staging" / "pubmed_gpt.txt"
REF_INHIBITORS = ROOT / "data" / "reference" / "mitochondrial_complex_i_inhibitors.txt"
BLACKLIST = ROOT / "data" / "reference" / "blacklist.txt"
PROCESSED_DIR = ROOT / "data" / "processed"

MCI_REFS = {
    "biguanide": "NC(=N)NC(=N)N",
    "metformin": "CN(C)C(=N)NC(=N)N",
    "phenformin": "N=C(N)NC(=N)NCCc1ccccc1",
    "buformin": "CCCCN=C(N)N=C(N)N", 
    "biguanide_motif": "N=C(N)NC(=N)NCCCCCCNC(=N)NC(=N)N",
    "proguanil": "CC(C)NC(=N)NC(=N)Nc1ccc(Cl)cc1"
}

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


def _make_morgan_gen(fp_radius: int, nbits: int):
    # Try modern kw signature first; fall back to positional for older builds.
    try:
        return rfg.GetMorganGenerator(
            radius=fp_radius,
            includeChirality=True,   # NOTE: includeChirality (not useChirality)
            useBondTypes=True,
            fpSize=nbits
        )
    except TypeError:
        # Older RDKit: use positional args to avoid kw-name drift
        # (radius, countSimulation, includeChirality, useBondTypes,
        #  onlyNonzeroInvariants, includeRingMembership, countBounds, fpSize, ...)
        return rfg.GetMorganGenerator(fp_radius, False, True, True, False, True, None, nbits)

def add_tanimoto_scores(df: pd.DataFrame, fp_radius=2, nbits=2048, topk=3) -> pd.DataFrame:
    gen = _make_morgan_gen(fp_radius, nbits)

    def mol_from_smiles(s):
        if not isinstance(s, str) or not s.strip():
            return None
        try:
            return Chem.MolFromSmiles(s, sanitize=True)
        except Exception:
            return None

    def ecfp(m):
        if gen is not None:
            return gen.GetFingerprint(m)
        # Fallback for very old RDKit (deprecated but safe as a backup)
        return AllChem.GetMorganFingerprintAsBitVect(m, radius=fp_radius, nBits=nbits)

    def inchikey(m):
        try:
            return inchi.MolToInchiKey(m)
        except Exception:
            return None

    def topk_mean(vals, k=3):
        if not vals:
            return np.nan
        arr = sorted(vals, reverse=True)
        return float(np.mean(arr[: min(k, len(arr))]))

    def confidence_from_sims(maxsim, topkmean):
        v = topkmean if (topkmean == topkmean) else maxsim
        if v != v or v is None:
            return None
        if v >= 0.70: return "high"
        if v >= 0.50: return "medium"
        if v >= 0.30: return "low"
        return "very-low"

    refs = df[(df["known_status"].astype(str).str.lower() == "known") & df["SMILES"].notna()].copy()
    refs["mol"] = refs["SMILES"].apply(mol_from_smiles)
    refs = refs.dropna(subset=["mol"]).copy()
    refs["inchikey"] = refs["mol"].apply(inchikey)
    refs = refs.drop_duplicates(subset=["inchikey"]).reset_index(drop=True)
    refs["fp"] = refs["mol"].apply(ecfp)

    if refs.empty:
        out = df.copy()
        out["MaxSim_all"] = None
        out["TopKMean_all"] = None
        out["BestRef_name"] = None
        out["confidence_similarity"] = None
        return out

    ref_fps   = refs["fp"].tolist()
    ref_names = refs["compound"].tolist()

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
    
# pip install rdkit-pypi pandas
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolops import RemoveHs

MCI_REFS = {
    "biguanide": "NC(=N)NC(=N)N",
    "metformin": "CN(C)C(=N)NC(=N)N",
    "phenformin": "N=C(N)NC(=N)NCCc1ccccc1",
    "buformin": "CCCCN=C(N)N=C(N)N",
    "biguanide_motif": "N=C(N)NC(=N)NCCCCCCNC(=N)NC(=N)N",
    "proguanil": "CC(C)NC(=N)NC(=N)Nc1ccc(Cl)cc1"
}

def score_biguanide_like(smiles_list, refs=MCI_REFS, alpha=0.7, beta=0.3):
    """
    Given a list of SMILES, return a DataFrame with:
      - biguanide/biguanide_motif substructure flags
      - similarity vs 'biguanide' (Tversky/Dice)
      - best similarity across all refs (Tversky/Dice) + which ref matched
      - ranked by substructure hit and Tversky score
    """
    def largest_fragment_smiles(smiles: str):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None
            frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
            if not frags:
                return None
            frag = max(frags, key=lambda m: m.GetNumHeavyAtoms())
            return Chem.MolToSmiles(frag, isomericSmiles=True)
        except Exception:
            return None

    def mol_from_smiles(smi: str):
        if not smi:
            return None
        smi2 = largest_fragment_smiles(smi)
        if not smi2:
            return None
        m = Chem.MolFromSmiles(smi2)
        if not m:
            return None
        try:
            Chem.SanitizeMol(m)
        except Exception:
            return None
        return RemoveHs(m)

    def ecfp4_count_fp(mol):
        return AllChem.GetMorganFingerprint(mol, radius=2)  # count-based fingerprint

    def tversky(fp_a, fp_b):
        return DataStructs.TverskySimilarity(fp_a, fp_b, alpha=alpha, beta=beta)

    def dice(fp_a, fp_b):
        return DataStructs.DiceSimilarity(fp_a, fp_b)

    # prepare reference molecules/fingerprints
    ref_mols = {name: mol_from_smiles(smi) for name, smi in refs.items()}
    ref_fps  = {name: ecfp4_count_fp(m) for name, m in ref_mols.items() if m is not None}

    patt_biguanide = Chem.MolFromSmarts(refs.get("biguanide")) if "biguanide" in refs else None
    patt_biguanide_motif = Chem.MolFromSmarts(refs.get("biguanide_motif")) if "biguanide_motif" in refs else None

    rows = []
    for smi in smiles_list:
        m = mol_from_smiles(smi)
        if m is None:
            rows.append({
                "smiles": smi,
                "has_biguanide_core": False,
                "has_biguanide_motif": False,
                "sim_biguanide_tversky": None,
                "sim_biguanide_dice": None,
                "best_biguanide_like_tversky": None,
                "best_ref_name_tversky": None,
                "best_biguanide_like_dice": None,
                "best_ref_name_dice": None
            })
            continue

        has_core = patt_biguanide and m.HasSubstructMatch(patt_biguanide)
        has_motif = patt_biguanide_motif and m.HasSubstructMatch(patt_biguanide_motif)

        fp = ecfp4_count_fp(m)

        sim_t = round(tversky(fp, ref_fps["biguanide"]), 3)
        sim_d = round(dice(fp, ref_fps["biguanide"]), 3)

        t_vals = {name: tversky(fp, rfp) for name, rfp in ref_fps.items()}
        d_vals = {name: dice(fp, rfp) for name, rfp in ref_fps.items()}
        best_t = max(t_vals, key=t_vals.get)
        best_d = max(d_vals, key=d_vals.get)

        rows.append({
            "smiles": smi,
            "has_biguanide_core": bool(has_core),
            "has_biguanide_motif": bool(has_motif),
            "sim_biguanide_tversky": sim_t,
            "sim_biguanide_dice": sim_d,
            "best_biguanide_like_tversky": round(t_vals[best_t], 3),
            "best_ref_name_tversky": best_t,
            "best_biguanide_like_dice": round(d_vals[best_d], 3),
            "best_ref_name_dice": best_d
        })

    df = pd.DataFrame(rows)
    df["_rank_key"] = df.apply(
        lambda r: (0 if r["has_biguanide_core"] or r["has_biguanide_motif"] else 1,
                   -(r["best_biguanide_like_tversky"] or -1)), axis=1)
    df = df.sort_values("_rank_key").drop(columns=["_rank_key"]).reset_index(drop=True)
    
    return df


# --- Example usage ---
# smiles_list = ["NC(=N)NC(=N)N", "CN(C)C(=N)NC(=N)N", "CCCN=C(N)N=C(N)N", "CC(C)NC(=N)NC(=N)Nc1ccc(Cl)cc1"]
# df = score_biguanide_like(smiles_list)
# print(df)

    
    
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
black_ref = set([e.strip().lower() for e in ref if e.strip()])

blacklist = []
with open(BLACKLIST, "r", encoding="utf-8") as f:
    blacklist = sorted(set([e.lower().strip() for e in f if e.strip()]))

blacklist = set(blacklist)
blacklist.update(black_ref)

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
ref.difference_update(["Roterone","Piericidin","Bongkrekic","IACS-10759"])
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
print("BLACKLSIT ->",sorted(blacklist))
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
    inh = [(e[0],e[1],e[2].replace("analogs","").replace("analogue","").replace("analog","").replace("diphenyleneiodonium","diphenylene iodonium").replace("acetogenins","acetogenin").replace("aroclor 1254","aroclor-1254").replace("annonaceous acetogenin","acetogenin").replace("deltalac-acetogenin","acetogenin")) for e in inh if e[2].strip()]
    inh = [(e[0],e[1],e[2].strip()) for e in inh if e[2].strip()]
    inh = [e for e in inh if e[1] and e[1].lower() != "no" and e[2] and e[2].lower() != "na"]
    inh = [e for e in inh if len(e[2]) > 2 and e[2].lower() not in blacklist]
    inh = [e for e in inh if e[2].lower()!="no" and e[2].lower().find("nitric oxide")==-1 and e[2].lower().find("mitochondr")==-1 and e[2].lower().find("silencing")==-1 and (not e[2].lower().startswith("compound")) and e[2].lower().find("nanoparticle")==-1 ]

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

stats["confidence_pubmed"] = pd.cut(
    stats["pubmed_references"],
    bins=[-np.inf, 1, 2, 4, np.inf],
    labels=["very-low", "low", "medium", "high"]
)

stats = stats[["compound", "pubmed_references", "known_status", "confidence_pubmed", "pubmed_ids"]]


# Final sort
stats = stats.sort_values(["pubmed_references", "compound"], ascending=[False, True]).reset_index(drop=True)

# add smiles
###############################################################################


#targets = stats.loc[stats["confidence_pubmed"]!="very-low","compound"].unique().tolist()
# take all compounds
targets = stats["compound"].unique().tolist()

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
        time.sleep(1)
        smi, source, query = fetch_smiles(t,normalize=False)          
    print(f"{i}/{n} compounds: {t} -> {smi} - {source} - {query}")  
    smiles[t] = smi
    time.sleep(1)

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
        "confidence_bins": "very-low: ≤1; low: 2; medium: 3–4; high: ≥5",
        "pubmed_ids_concat": "unique PMIDs, sorted, ';' separated",
    },
    notes=(
        "Aggregated stats per compound including:\n"
        "- pubmed_references (unique PMID count)\n"
        "- known_status (reference list appended at 100)\n"
        "- confidence_pubmed: categorical score based on pubmed_references count "
        "(very-low, low, medium, high)\n"
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
        "confidence_bins": "very-low: ≤1; low: 2; medium: 3–4; high: ≥5",
        "pubmed_ids_concat": "unique PMIDs, sorted, ';' separated",
    },
    notes=(
        "Excel mirror of the aggregated compound summary with SMILES and RDKit "
        "similarity metrics."
    ),
)

print('[INFO] Setup complete.')
