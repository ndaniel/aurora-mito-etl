#!/usr/bin/env bash
# Aurora ETL data downloader
# Automatically detects its path and downloads PubTator + PubMed data
# into ../data/raw relative to where this script is located.

set -euo pipefail
IFS=$'\n\t'

# ---------------------------------------------------------------------------
# Detect paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$(cd "$SCRIPT_DIR/../data/raw" && pwd)"

cd "$DATA_DIR"
echo "Script location: $SCRIPT_DIR"
echo "Data directory:  $DATA_DIR"

# download PubTator
echo "Downloading PubTator"
mkdir -p pubtator
wget https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/chemical2pubtator3.gz -O pubtator/chemical2pubtator3.gz
curl -s https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/ | grep chemical2pubtator3.gz | awk '{print $3}' > pubtator/release_info.txt
curl -s https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/ | grep chemical2pubtator3.gz >> pubtator/release_info.txt

# download PubMed
echo "Downloading Pubmed - abstracts"
mkdir -p pubmed
curl -s https://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/ | grep .xml.gz | grep -v .xml.gz.md5 | head -1 | awk '{print $3}' > pubmed/release_info.txt
wget -nd -r -l1 -np -A xml.gz -c -N -e robots=off --directory-prefix="pubmed" https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/
wget -nd -r -l1 -np -A xml.gz -c -N -e robots=off --directory-prefix="pubmed" https://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/

# Download a list of MESH terms that could be used as inhibitor (small compounds, organic, etc.)
echo "Downloading Pubmed - MESH info"
mkdir -p mesh
rm -f mesh/desc2025.xml
rm -f mesh/supp2025.xml
wget --directory-prefix="mesh" https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2025.xml 
wget --directory-prefix="mesh" https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/supp2025.xml 









