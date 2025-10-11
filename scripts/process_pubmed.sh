#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

# ------------------------------------------------------------------------------
# PubMed XML → TSV: PMID, YEAR, (TITLE + '.' if needed + ' ' + ABSTRACT)
# - One UTF-8→ASCII-cleaned output file per input XML.GZ
# - Steps:
#     1. Extract with xmlstarlet
#     2. AWK filters (inhibitor + mito-complex OR NADH-enzyme)
#     3. Strip inline tags
#     4. Unicode cleanup (NFC + remove NBSP, zero-widths, BOM)
#     5. Drop control chars
#     6. Transliterate to ASCII (iconv //TRANSLIT)
# - Fails if uconv or iconv is missing
# ------------------------------------------------------------------------------

# Required tools
command -v xmlstarlet >/dev/null
command -v parallel   >/dev/null
command -v gzip       >/dev/null
command -v uconv      >/dev/null
command -v iconv      >/dev/null

# Threads
NP="$( (command -v nproc >/dev/null && nproc) || echo 2 )"
DEF_THREADS=$(( NP > 1 ? NP - 1 : 1 ))
THREADS="${1:-$DEF_THREADS}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$(cd "$SCRIPT_DIR/../data/staging" && pwd)"
INROOT="$DATA_DIR/../raw/pubmed"
OUTDIR="$DATA_DIR/pubmed"
BIGFILE="$DATA_DIR/pubmed.txt"
YEAR_MIN=2000
mkdir -p "$OUTDIR"

echo "Input root: $INROOT"
echo "Output dir: $OUTDIR"
echo "Threads:    $THREADS"

find "$OUTDIR" -maxdepth 1 -type f -name 'pubmed*.txt' -delete 2>/dev/null || true
find "$OUTDIR" -maxdepth 1 -type f -name '.*.tmp*'    -delete 2>/dev/null || true

DECOMP="${DECOMP:-gzip -dc}"
TAB="$(printf '\t')"
export OUTDIR DECOMP TAB

# Year XPath
DATE_XPATH='normalize-space((
  MedlineCitation/DateCompleted/Year |
  PubmedData/History/PubMedPubDate[@PubStatus="pubmed"]/Year |
  MedlineCitation/DateRevised/Year
)[1])'
export DATE_XPATH

# ------------------------------------------------------------------------------
# Regex filters (POSIX ERE; awk-friendly)
# ------------------------------------------------------------------------------
RX_INHIBIT='(inhibit|antagoni|block|down[ -]?regulat|impair|repress)'

# "mitochondria ... complex" OR "complex ... mitochondria" (order-agnostic, lightweight)
RX_COMPLEX='(mitochond[^[:space:]]*.*complex|complex[^[:space:]]*.*mitochond)'

# NADH … (dehydrogenase|oxidoreductase), either order, allowing punctuation/words in between
RX_NADH='(NADH([[:space:][:punct:]]+[^[:space:]]+){0,8}[[:space:][:punct:]]+(dehydrogenase|oxidoreductase)|(dehydrogenase|oxidoreductase)([[:space:][:punct:]]+[^[:space:]]+){0,8}[[:space:][:punct:]]+NADH)'

export RX_INHIBIT RX_COMPLEX RX_NADH

find "$INROOT" -type f -name '*.xml.gz' -print0 | sort -z |
parallel -0 -j"$THREADS" --halt now,fail=1 --linebuffer '
  infile="{}"
  base="$(basename "$infile" .xml.gz)"
  outfile="$OUTDIR/${base}.txt"
  tmp_raw="$OUTDIR/.${base}.tmp.raw"
  tmp_mid="$OUTDIR/.${base}.tmp.mid"
  tmp_nfc="$OUTDIR/.${base}.tmp.nfc"
  tmp_utf="$OUTDIR/.${base}.tmp.utf"
  tmpfile="$OUTDIR/.${base}.tmp"

  echo "[INFO] Processing $infile → $outfile"

  # XML extract
  '"$DECOMP"' "$infile" | tail -n +3 |
  xmlstarlet sel --nonet --nocdata -T -t \
    -m "//PubmedArticle" \
    -v "(MedlineCitation/PMID[@Version=\"1\"] | MedlineCitation/PMID)[1]" -o "$TAB" \
    -v "$DATE_XPATH" -o "$TAB" \
    -v "normalize-space(string(MedlineCitation/Article/ArticleTitle))" -o "$TAB" \
    -v "normalize-space(string(MedlineCitation/Article/Abstract))" \
    -n 2>/dev/null |
  awk -F'\''\t'\'' -v OFS='\''\t'\'' -v IGNORECASE=1 \
      -v RXI="$RX_INHIBIT" -v RXC="$RX_COMPLEX" -v RXN="$RX_NADH" "
  {
    for (i=1; i<=NF; i++) {
      gsub(/[[:space:]]+/, \" \", \$i)
      sub(/^[[:space:]]+/, \"\", \$i)
      sub(/[[:space:]]+$/, \"\", \$i)
    }
    pmid=\$1; raw_year=\$2; title=\$3; abstract=\$4
    if (title==\"\" || abstract==\"\") next
    if (title ~ /^\[/) next
    if (!match(raw_year, /(19|20)[0-9]{2}/)) next
    yr=substr(raw_year, RSTART, 4)
    t=title; if (t !~ /\.$/) t=t\".\"
    ta=t\" \"abstract

    has_inhib      = (ta ~ RXI)
    is_mito_complex= (ta ~ RXC)
    is_nadh_enzyme = (ta ~ RXN)

    if (has_inhib && (is_mito_complex || is_nadh_enzyme))
      print pmid, yr, ta
  }" > "$tmp_raw"

  # Remove inline tags
  sed -E "s#</?(sup|sub|i|b|u|em|strong|small|p|br)\\b[^>]*>##gI" "$tmp_raw" > "$tmp_mid"

  # Normalize Unicode (NFC)
  uconv -x "Any-NFC" < "$tmp_mid" > "$tmp_nfc"

  # Replace NBSP, remove zero-width & BOM
  NBSP="$(printf '\''\302\240'\'')"     # U+00A0
  ZWSP="$(printf '\''\342\200\213'\'')" # U+200B
  ZWNJ="$(printf '\''\342\200\214'\'')" # U+200C
  ZWJ="$(printf '\''\342\200\215'\'')"  # U+200D
  BOM="$(printf '\''\357\273\277'\'')"  # U+FEFF
  sed "s/${NBSP}/ /g; s/${ZWSP}//g; s/${ZWNJ}//g; s/${ZWJ}//g; s/${BOM}//g" \
      "$tmp_nfc" > "$tmp_utf"

  # Drop control chars
  tr -d '\''\000-\010\013\014\016-\037\177'\'' < "$tmp_utf" > "$tmpfile"

  # Transliterate to ASCII (mandatory)
  iconv -f UTF-8 -t ASCII//TRANSLIT < "$tmpfile" > "$outfile"

  # Cleanup
  rm -f "$tmp_raw" "$tmp_mid" "$tmp_nfc" "$tmp_utf" "$tmpfile"

  echo "[DONE]  $infile → $outfile"
'

# --------------------------------------------------------------------------
# Concatenate all cleaned PubMed TXT files (ascending by name)
# Keep only lines where YEAR (2nd column) >= YEAR_MIN
# --------------------------------------------------------------------------
rm -f "$BIGFILE"

find "$OUTDIR" -maxdepth 1 -type f -name 'pubmed*.txt' -print0 | sort -z | \
xargs -0 -r cat -- | \
awk -F'\t' -v OFS='\t' -v Y="$YEAR_MIN" '$2 ~ /^[0-9]{4}$/ && ($2+0) >= Y' > "$BIGFILE"

LINES=$(wc -l < "$BIGFILE" | tr -d " ")
echo "[DONE]  Final PubMed file contains $LINES lines → $BIGFILE"

# Remove the intermediate directory with per-file outputs
echo "[INFO] Cleaning up intermediate files: $OUTDIR"
rm -rf "$OUTDIR"
echo "[DONE]  Removed $OUTDIR"




