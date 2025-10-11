#!/usr/bin/env bash

# Run the processing for staging area files

./process_mesh.py
./process_pubtator.sh
./process_pubmed.sh
./merge_filter.sh
./run_gpt_filter.py





