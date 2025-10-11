#!/usr/bin/env python3
"""
Build a broad list of bioactive small-molecule MeSH concepts (Descriptors + SCRs).

- Descriptors (D-IDs): keep if any TreeNumber starts with D02–D06 (organic-ish).
- Supplemental Concept Records (C-IDs): keep if
    (a) maps to a kept Descriptor (D02–D06) OR to a pharmacologic-action branch (D27.*), OR
    (b) has a CAS-like RegistryNumber (NNNNN-NN-N), OR
    (c) is an SCRClass="1" pharmacologic substance with a code-name (e.g., IACS-010759, ODM-208) and RN ∈ {"", "0"}.

Outputs:
    data/staging/mesh_bioactive.tsv         (Type, MeSH_UI, Name, OneTreeNumber)
    data/staging/mesh_bioactive_tags.txt    (MESH:D..., MESH:C...)
"""

import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

# --------- Paths (repo-layout aware) ----------
ROOT = Path(__file__).resolve().parents[1]  # repo root (etl/… → repo/)
DESC_XML = ROOT / "data" / "raw" / "mesh" / "desc2025.xml"
SUPP_XML = ROOT / "data" / "raw" / "mesh" / "supp2025.xml"
OUT_TSV  = ROOT / "data" / "staging" / "mesh_bioactive.tsv"
OUT_TAGS = ROOT / "data" / "staging" / "mesh_bioactive_tags.txt"

# --------- Config ----------
# Organic/small-molecule–rich branches (Descriptors). Widen if needed (e.g., include D08 for lipids).
ORGANIC_BRANCH = re.compile(r"^D0(2|3|4|5|6)\.")
# Pharmacologic action branch (Descriptors) used to widen SCR mapping only.
PHARMA_BRANCH = re.compile(r"^D27\.")
# Simple CAS number pattern (NNNNN-NN-N)
CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")
# Strip leading '*' or spaces from mapped DescriptorUI (major topic marker).
CLEAN_DUI = re.compile(r'^[*\s]+')

# Heuristic for code-named small molecules
# Examples matched: IACS-010759, GDC-0941, BAY 11-7082, AZD9291, ODM-208, ORM-12741, MK2206, BMS986205, ONO-7475
CODE_NAME_RE = re.compile(
    r"""(?ix)
    ^(
        [A-Z]{2,5}[A-Z0-9]*          # sponsor/code prefix (2–5 uppercase letters, may include digits)
        (?:[-\s]?\d{2,6}[A-Z0-9]*)+  # one or more numeric chunks (allows hyphen/space)
    )$
    """
)
# ---------------------------

def t(elem) -> str:
    return (elem.text or "").strip() if elem is not None else ""

def require_file(path: Path, download_url: str | None = None):
    if path.exists() and path.stat().st_size > 0:
        return
    if not download_url:
        sys.stderr.write(f"Missing {path}\n")
        sys.exit(1)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import urllib.request
        urllib.request.urlretrieve(download_url, path)
        if not path.exists() or path.stat().st_size == 0:
            raise RuntimeError("download produced empty file")
    except Exception as e:
        sys.stderr.write(f"Failed to fetch {download_url} → {path} : {e}\n")
        sys.exit(1)

def read_descriptors(path: Path):
    """
    Returns dict:
      dui -> {"name": str, "trees": [str]}
    """
    d = {}
    root = ET.parse(path).getroot()
    for rec in root.findall("DescriptorRecord"):
        dui = t(rec.find("DescriptorUI"))
        if not dui:
            continue
        name = t(rec.find("./DescriptorName/String"))
        trees = [t(x) for x in rec.findall("./TreeNumberList/TreeNumber")]
        d[dui] = {"name": name, "trees": trees}
    return d

def read_scrs(path: Path):
    """
    Returns:
      scr_names: CUI -> preferred name (SupplementalRecordName/String)
      scr_cas:   CUI -> RegistryNumber (may be '0' or '')
      scr_map:   CUI -> set of mapped Descriptor UIs (major-topic '*' removed)
      scr_class: CUI -> SCRClass attribute (e.g., '1' for Pharmacologic Substance)
    """
    scr_names, scr_cas, scr_map, scr_class = {}, {}, {}, {}
    root = ET.parse(path).getroot()
    for rec in root.findall("SupplementalRecord"):
        cui = t(rec.find("SupplementalRecordUI"))
        if not cui:
            continue
        scr_names[cui] = t(rec.find("./SupplementalRecordName/String"))
        scr_cas[cui]   = t(rec.find("./RegistryNumber"))
        scr_class[cui] = rec.get("SCRClass", "")

        mapped: list[str] = []
        for h in rec.findall("./HeadingMappedToList/HeadingMappedTo"):
            dui = t(h.find("./DescriptorReferredTo/DescriptorUI")) or t(h.find("./DescriptorUI"))
            if dui:
                dui = CLEAN_DUI.sub("", dui)  # remove leading '*' / spaces
                mapped.append(dui)
        if mapped:
            scr_map[cui] = set(mapped)
    return scr_names, scr_cas, scr_map, scr_class

def descriptor_is_organic(desc_meta) -> bool:
    return any(ORGANIC_BRANCH.search(tr) for tr in desc_meta["trees"])

def descriptor_is_pharma_or_organic(desc_meta) -> bool:
    return any(ORGANIC_BRANCH.search(tr) or PHARMA_BRANCH.search(tr) for tr in desc_meta["trees"])

def main():
    # Ensure inputs exist (auto-download if missing)
    require_file(DESC_XML, "https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2025.xml")
    require_file(SUPP_XML, "https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/supp2025.xml")

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)

    desc = read_descriptors(DESC_XML)
    scr_names, scr_cas, scr_map, scr_class = read_scrs(SUPP_XML)

    # Kept descriptors for output (strict organic branches)
    keepD: dict[str, bool] = {dui: descriptor_is_organic(meta) for dui, meta in desc.items()}

    # Kept descriptors for SCR mapping (organic OR pharmacologic-action)
    keep_for_mapping: dict[str, bool] = {dui: descriptor_is_pharma_or_organic(meta) for dui, meta in desc.items()}

    rows: list[tuple[str, str, str, str]] = []

    # 1) Add qualifying Descriptors (D-IDs)
    for dui, meta in desc.items():
        if keepD.get(dui):
            one_tree = meta["trees"][0] if meta["trees"] else ""
            rows.append(("D", dui, meta["name"], one_tree))

    # 2) Add SCRs (C-IDs) per widened logic
    for cui, name in scr_names.items():
        mapped_ds = scr_map.get(cui, set())
        maps_to_kept_or_pharma = any(keep_for_mapping.get(d, False) for d in mapped_ds)
        rn = scr_cas.get(cui, "")
        has_cas = bool(CAS_RE.match(rn))
        looks_like_code_name = bool(CODE_NAME_RE.match(name))
        is_pharm_substance = scr_class.get(cui, "") == "1"

        if maps_to_kept_or_pharma or has_cas or (is_pharm_substance and looks_like_code_name and rn in ("", "0")):
            # pick one tree to display (from any mapped kept-or-pharma descriptor)
            one_tree = ""
            for d in mapped_ds:
                if keep_for_mapping.get(d, False):
                    trees = desc.get(d, {}).get("trees", [])
                    if trees:
                        one_tree = trees[0]
                        break
            rows.append(("C", cui, name, one_tree))

    # 3) De-duplicate and write outputs
    seen = set()
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w", encoding="utf-8") as f:
        f.write("Type\tMeSH_UI\tName\tOneTreeNumber\n")
        for typ, ui, name, one_tree in rows:
            key = (typ, ui)
            if key in seen:
                continue
            seen.add(key)
            f.write(f"{typ}\t{ui}\t{name}\t{one_tree}\n")

    # One tag per unique UI
    tag_set = {f"MESH:{ui}\n" for (_typ, ui) in {(typ, ui) for typ, ui, *_ in rows}}
    tags = sorted(set(tag_set))
    with OUT_TAGS.open("w", encoding="utf-8") as f:
        f.writelines(tags)

    print(f"Wrote {len(seen)} rows -> {OUT_TSV}")
    print(f"Tags -> {OUT_TAGS}")
    print("Sanity checks:")
    print(f"  grep -P '\\tC071942\\t' {OUT_TSV}   # arctigenin by ID")
    print(f"  grep -i 'arctigenin' {OUT_TSV}      # arctigenin by name (if preferred name matches)")
    print(f"  grep -P '\\tIACS-010759\\t' {OUT_TSV} # code-named example")
    print("Tip: widen ORGANIC_BRANCH (e.g., include D08) if needed.")
    print("     You can later intersect OUT_TAGS with PubTator column 3 to get PMIDs.")

if __name__ == "__main__":
    main()

