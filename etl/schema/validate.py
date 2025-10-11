#!/usr/bin/env python3
"""
Validate TSV/CSV artifacts against JSON Schemas in etl/schema/*.schema.json.

Usage examples:
  python etl/schema/validate.py data/staging/pubmed_gpt.txt etl/schema/pubmed_gpt.schema.json --tsv
  python etl/schema/validate.py data/processed/2025-10-11/all_mito_complex_I_inhibitors.txt etl/schema/processed_all.schema.json --tsv
"""
import argparse, json, sys
from pathlib import Path

def load_schema(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))

def parse_table(path: Path, is_tsv: bool):
    import csv
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t" if is_tsv else ",")
        for row in reader:
            yield row

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table", type=Path)
    ap.add_argument("schema", type=Path)
    ap.add_argument("--tsv", action="store_true", help="Input is TSV (default CSV)")
    args = ap.parse_args()

    try:
        import jsonschema
    except ImportError:
        print("[ERROR] jsonschema not installed. Run: pip install jsonschema", file=sys.stderr)
        sys.exit(2)

    schema = load_schema(args.schema)
    validator = jsonschema.Draft202012Validator(schema)

    errors = 0
    for i, row in enumerate(parse_table(args.table, is_tsv=args.tsv), start=1):
        for err in validator.iter_errors(row):
            print(f"[ERROR] row {i}: {err.message}", file=sys.stderr)
            errors += 1

    if errors:
        print(f"[FAIL] {errors} validation error(s).", file=sys.stderr)
        sys.exit(1)
    else:
        print("[OK] Validation passed.")

if __name__ == "__main__":
    main()
