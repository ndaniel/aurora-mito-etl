#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

# ------------------------------------------------------------------------------
# Step 1: Remove all lines from pubmed.txt that contain any inhibitor terms
#          (case-insensitive, fast with ripgrep)
# Step 2: Keep only lines whose PubMed ID (col 1) exists in pubtator_pmids.txt
#
# Multithreaded:
#   - Detects CPU cores and uses (nproc - 1)
#   - sort uses --parallel and temporary files on TMPDIR
#
# Input:
#   data/staging/pubmed.txt
#   data/reference/mitochondrial_complex_i_inhibitors.txt
#   data/staging/pubtator_pmids.txt
#
# Output:
#   data/staging/pubmed_filtered.txt
# ------------------------------------------------------------------------------

# Detect hardware threads (at least 1)
NP="$( (command -v nproc >/dev/null && nproc) || echo 2 )"
THREADS=$(( NP > 1 ? NP - 1 : 1 ))

DATA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../data" && pwd)"
STAGING="$DATA_DIR/staging"
REF="$DATA_DIR/reference"

INPUT="$STAGING/pubmed.txt"
INHIBITORS_RAW="$REF/mitochondrial_complex_i_inhibitors.txt"
PUBTATOR_PMIDS="$STAGING/pubtator_pmids.txt"
TMP1="$STAGING/.pubmed_tmp1.txt"
TMP2="$STAGING/.pubtator_pmids_sorted.txt"
TMP3="$STAGING/.pubmed_filtered.txt"
TMP_KEYED="$STAGING/.pubmed_bykey.txt"
INHIBITORS_CLEAN="$STAGING/.inhibitors.clean.txt"
OUTPUT="$STAGING/pubmed_filtered.txt"

echo "[INFO] Detected $NP cores → using $THREADS threads for sorting"
echo "[INFO] Input: $INPUT"
echo "[INFO] Inhibitors: $INHIBITORS_RAW"
echo "[INFO] PubTator PMIDs: $PUBTATOR_PMIDS"
echo "[INFO] Output: $OUTPUT"

# ---------------------------- Checks -----------------------------------------
for f in "$INPUT" "$INHIBITORS_RAW" "$PUBTATOR_PMIDS"; do
  [[ -f "$f" ]] || { echo "[ERROR] Missing: $f" >&2; exit 1; }
done
command -v rg >/dev/null || { echo "[ERROR] ripgrep (rg) not found" >&2; exit 1; }
command -v sort >/dev/null || { echo "[ERROR] sort not found" >&2; exit 1; }
command -v join >/dev/null || { echo "[ERROR] join not found" >&2; exit 1; }

# --------------------------- Step 0: clean patterns ---------------------------
# Strip comments/blanks and very-short tokens (length<3) that overmatch.
sed -e 's/#.*$//' -e 's/^[[:space:]]\+//' -e 's/[[:space:]]\+$//' \
    -e '/^$/d' "$INHIBITORS_RAW" | awk 'length>=3' > "$INHIBITORS_CLEAN"

# --------------------------- Step 1 ------------------------------------------
echo "[INFO] Step 1: Removing lines matching inhibitors (fixed strings, -f)…"
rg -v -i -F -f "$INHIBITORS_CLEAN" "$INPUT" > "$TMP1"

LINES1=$(wc -l < "$TMP1" | tr -d " ")
echo "[INFO]  → Remaining after inhibitor filtering: $LINES1 lines"

# --------------------------- Step 2 ------------------------------------------
echo "[INFO] Step 2: Keeping only lines whose PubMed ID (col1) is in pubtator_pmids.txt"

# Sort both inputs (TMP1 + TMP2)
sort --parallel="$THREADS" -t$'\t' -k1,1 "$TMP1" -o "$TMP1"
sort --parallel="$THREADS" -u "$PUBTATOR_PMIDS" -o "$TMP2"

# Prepare TMP1 for join: prefix with PMID as key (for join field 1)
awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$0}' "$TMP1" > "$TMP_KEYED"
sort --parallel="$THREADS" -t$'\t' -k1,1 "$TMP_KEYED" -o "$TMP_KEYED"

# Join intersection on PMID (field 1)
join -t $'\t' -1 1 -2 1 "$TMP_KEYED" "$TMP2" \
  | cut -f2- \
  | sort --parallel="$THREADS" -t$'\t' -k1,1 -o  "$TMP3"

rg -iP '(?:(?=.*\bmitochondr)\bcomplex(?:es)?\b[\s\-]*(?:(?:I|1)(?:\s*[-–—−]\s*(?:II|2|III|3|IV|4|V|5))?)\b|\bNADH\b(?:\s*[:\-–—−]\s*|[-\s]*(?:dependent|linked)\s+)?(?:dehydrogenase(?:s)?|(?:ubiquinone|quinone)\s*oxidoreductase)\b)' \
  "$TMP3" > "$OUTPUT"

LINES2=$(wc -l < "$OUTPUT" | tr -d " ")
echo "[DONE]  Final filtered output → $OUTPUT"
echo "[INFO]  Remaining lines after PubTator PMID filter: $LINES2"

# Cleanup
rm -f "$TMP1" "$TMP2" "$TMP3" "$TMP_KEYED" "$INHIBITORS_CLEAN"
echo "[INFO] Cleanup complete."
echo "[INFO] All done."

