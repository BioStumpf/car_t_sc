#!/bin/bash

INPT_DIR="$HOME/car_t_sc/01_data/processed/OVA_receptor/hmmersearch_alignment_outputs"
OUTPT_DIR="$HOME/car_t_sc/01_data/processed/OVA_receptor/hmmersearch_alignment_outputs/ids_txts_only"
SCRIPT="$HOME/car_t_sc/02_code/OVA_mapping/bowtie2_of_hmmer_hits/extract_ids_to_txt.py"

for file in $INPT_DIR/*.out; do
    pool_pattern=$(basename $file | grep -o "P[1-9]")
    lib_pattern=$(basename $file | grep -oE "GEX|VDJ")
    output_file="${OUTPT_DIR}/${pool_pattern}_${lib_pattern}.txt"
    python3 "$SCRIPT" -hmmer "$file" -out "$output_file"
done