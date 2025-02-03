#!/bin/bash

ANNOTATED_DIR="$HOME/car_t_sc/01_data/processed/OVA_receptor/hmmer_hits_bowtie2/annotated_genes_for_bowtie_hits"
HMMER_HITS_DIR="$HOME/car_t_sc/01_data/processed/OVA_receptor/hmmersearch_alignment_outputs"
OUTPT_DIR="$HOME/car_t_sc/01_data/processed/OVA_receptor/hmmersearch_alignment_outputs/non_tcr_filtered_tabular_outs"
SCRIPT="$HOME/car_t_sc/02_code/OVA_mapping/filtering_hmmer_hits/filtering.py"

for file in $ANNOTATED_DIR/*.bed; do
    pool_pattern=$(basename $file | grep -o "P[1-9]")
    lib_pattern=$(basename $file | grep -oE "GEX|VDJ")
    output_file="${OUTPT_DIR}/${pool_pattern}_${lib_pattern}.out"
    hmmer_hits="${HMMER_HITS_DIR}/${pool_pattern}_${lib_pattern}.out"
    python3 $SCRIPT -hmmer_annotated $file -hmmer $hmmer_hits -out $output_file
done