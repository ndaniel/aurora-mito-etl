# ETL Data Dictionary

This document defines the canonical schemas for the Aurora MCI ETL pipeline.
Column order is normative for TSV/CSV artifacts.

## Staging

### `data/staging/pubmed_gpt.txt`
- **pmid** *(string)* — PubMed ID (digits).
- **confidence** *(string)* — one of `YES`, `probablyYES`, `NO`.
- **compound** *(string)* — small molecule name (free text).

### `data/staging/mesh_bioactive.tsv`
- **Type** *(string)* — `Descriptor` or `SCR`.
- **MeSH_UI** *(string)* — MeSH ID (e.g., `D012345`, `C123456`).
- **Name** *(string)* — preferred concept name.
- **OneTreeNumber** *(string)* — any one representative TreeNumber kept for the row.

### `data/staging/pubtator_filtered.tsv`
- **pmid** *(string)* — PubMed ID.
- **mention** *(string)* — surface form of the chemical.
- **normalized_id** *(string)* — normalized identifier (e.g., MeSH UI, CAS, or other).

## Processed

### `data/processed/<date>/new_mito_complex_I_inhibitors.txt`
- **pmid** *(string)* — PubMed ID.
- **confidence** *(string)* — `YES` or `probablyYES` (NOs omitted).
- **compound** *(string)* — predicted inhibitor mention.

### `data/processed/<date>/all_mito_complex_I_inhibitors.txt`
- **compound** *(string)* — canonical representative string.
- **pubmed_references** *(integer)* — count of unique PMIDs referencing this compound.
- **known_status** *(string)* — `known` (from reference list) or `new` (from classifier outputs).

## Conventions
- All files are UTF-8 encoded.
- TSV delimiter is `\t`. CSV not used unless explicitly stated.
- Header row is present in every artifact.
- `release_info.txt` sits next to outputs and lists SHA256, date, step, and sources.
