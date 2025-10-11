#!/usr/bin/env bash

set -euo pipefail
export LC_ALL=C
IFS=$'\n\t'

# ---------------------------------------------------------------------------
# Detect paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$(cd "$SCRIPT_DIR/../data/staging" && pwd)"

cd "$DATA_DIR"
echo "Script location: $SCRIPT_DIR"
echo "Data directory:  $DATA_DIR"

#awk -F '\t' '{ if ($3 != "") print $1 }' chemical2pubtator3 | sort -u> pubtator_pmids.txt

#zcat ../raw/pubtator/chemical2pubtator3.gz | awk -F '\t' '$3!=""{print $1}' | sort -n -u > pubtator_pmids.txt

# zcat ../raw/pubtator/chemical2pubtator3.gz \
#   | awk -F'\t' 'NR==FNR{ok[$2]=1; next} ($3 in ok){print $1}' \
#       mesh_bioactive.tsv - \
#   | sort -n -u > pubtator_pmids.txt

#sort -u data/staging/mesh_bioactive_tags.txt -o data/staging/mesh_bioactive_tags_sorted.txt

# 2. Join PubTator with the tag list and extract PMIDs
zcat ../raw/pubtator/chemical2pubtator3.gz \
  | awk -F'\t' '$3 != "" {print $3 "\t" $1}' \
  | sort -k1,1 -S1G \
  | join -t $'\t' -j 1 - mesh_bioactive_tags.txt \
  | cut -f2 \
  | sort -u > pubtator_pmids.txt

echo "Extracted $(wc -l < pubtator_pmids.txt) unique PubMed IDs from PubTator"
echo "Sample PubMed IDs:"
head -5 pubtator_pmids.txt  
echo "..."
tail -5 pubtator_pmids.txt
echo ""
echo "Done"
