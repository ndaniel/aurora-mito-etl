# AURORA-MITO ETL pipeline


A reproducible ETL (Extract–Transform–Load) pipeline for building a local,
queryable mirror of [PubMed](http://pubmed.ncbi.nlm.nih.gov/) abstracts 
and [PubTator](https://www.ncbi.nlm.nih.gov/research/pubtator3/) chemical annotations.

This project extracts small compounds names that are known to inhibit mitochondrial 
complex I from PubMed database using LLM and enriches them with cheminformatics
metadata for downstream analysis.

This repository automates:
1. **Data acquisition** – downloads the latest PubMed XML baselines + updates
   and PubTator3 chemical annotations.
2. **Transformation** – parses, normalizes, and converts the XML / TSV data
   into tables.
3. **Loading / Curation** – prepares analysis-ready datasets for downstream
   machine-learning or knowledge-graph projects.


---

## Quick start

### 1. Clone the repository

```bash
git clone https://github.com/ndaniel/aurora-mito-etl.git
cd aurora-mito-etl
```

### 2. Download data

```bash
bash scripts/run_pipeline.sh
```

This downloads and verifies:

- PubTator: `chemical2pubtator3.gz`
- PubMed: baseline and updatefiles (`*.xml.gz`)

The results are stored under `data/raw/`.

---

## System requirements

This pipeline requires a **GNU/Linux environment** (tested on Ubuntu 24.04 LTS).

Essential command-line tools:
- `bash` ≥ 4.0  
- `wget`
- `curl`
- `awk` (`gawk`)
- `parallel`
- `xmlstarlet`
- `pigz`
- `iconv`
- `uconv`
- `gzip`
- `ripgrep`
- `grep` (GNU)
- `zcat`

Python libraries are managed via `requirements.txt`. Recent updates to the
release finalization step require the following additional packages:
- `requests` – REST calls to PubChem and ChEMBL.
- `rdkit` – molecular fingerprints + Tanimoto similarity scoring.
- `openpyxl` – Excel writer backend for the per-run summary workbook.

---

## Provenance

Each data source stores its own metadata file:

- `data/raw/pubtator/release_info.txt`
- `data/raw/pubmed/release_info.txt`

These include:
- source URLs
- fetch timestamps
- file counts
- last-modified headers
- checksums

Ensuring every ETL run is fully traceable and reproducible.


---

## Data Sources and Licensing

This project **does not host, redistribute, or modify** PubMed or PubTator content.
All data are **downloaded directly** from official NCBI FTP servers:

- **PubMed**: [https://ftp.ncbi.nlm.nih.gov/pubmed/](https://ftp.ncbi.nlm.nih.gov/pubmed/)  
- **PubTator3**: [https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/](https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/)
- **PubChem**: [https://pubchem.ncbi.nlm.nih.gov/](https://pubchem.ncbi.nlm.nih.gov/)
- **ChEMBL**: [https://www.ebi.ac.uk/chembl/](https://www.ebi.ac.uk/chembl/)

The use, redistribution, and citation of the data are governed by their respective
providers - primarily the **U.S. National Library of Medicine (NLM)** 
and **National Center for Biotechnology Information (NCBI)**.

For details on licensing and reuse, refer to:

- [PubMed Data Usage Policies](https://www.nlm.nih.gov/databases/download/pubmed_medline.html)  
- [PubTator3 Terms of Use](https://www.ncbi.nlm.nih.gov/home/about/policies/)
- [PubChem Terms of Use](https://www.ncbi.nlm.nih.gov/home/about/policies/)
- [ChEMBL Terms of Use](https://chembl.github.io/chembl-licensing/)

All users of this pipeline are responsible for complying with the terms and
conditions of those data providers.

---

## Processed outputs

The finalization script (`scripts/finalize_realease.py`) assembles the
release artifacts under `data/processed/<date>/` and augments the
compound summary with cheminformatics diagnostics:

- Pulls SMILES strings via PubChem first, then falls back to ChEMBL.
- Generates RDKit Morgan fingerprints (ECFP4/2048) for similarity scoring.
- Assigns PubMed-driven confidence bins and complementary RDKit similarity labels.
- Flags biguanide-like chemistry to spotlight core pharmacophore analogs.
- Exports TSV artifacts for pipelines and Excel mirrors for analyst review.

The primary table `all_mito_complex_I_inhibitors.txt` now includes:
- Core attributes: `compound`, `pubmed_references`, `known_status`, `confidence_pubmed`, `pubmed_ids`.
- Similarity scores: `MaxSim_all`, `TopKMean_all`, `BestRef_name`, `confidence_similarity`.
- Biguanide diagnostics: `has_biguanide_core`, `has_biguanide_motif`, `sim_biguanide_tversky`, `sim_biguanide_dice`, `best_biguanide_like_tversky`, `best_ref_name_tversky`, `best_biguanide_like_dice`, `best_ref_name_dice`.
- Structural context: trailing `SMILES` column for downstream QSAR work.

An Excel mirror (`all_mito_complex_I_inhibitors.xlsx`) is emitted alongside the
TSV to simplify exploratory review.

Refer to `etl/schema/DATA_DICTIONARY.md` for the full column reference.


