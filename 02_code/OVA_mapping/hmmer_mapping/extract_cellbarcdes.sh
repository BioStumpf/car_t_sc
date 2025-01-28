#!/bin/bash

#specify Data and output directories
DATA_DIR="$HOME/car_t_sc/01_data"
FASTQS_DIR="$DATA_DIR/raw/raw/2024-06-17_Maik_Luu_24054SC_raw_FASTQ/FASTQ_by_lib"
OUTPUT="$DATA_DIR/processed/OVA_receptor/OVA_counts"
SCRIPT="$HOME/car_t_sc/02_code/OVA_mapping/hmmer_mapping/extract_cellbarcodes.py"
ALIGNMENT_DIR="$DATA_DIR/processed/OVA_receptor/hmmersearch_alignment_outputs"

for file in $FASTQS_DIR/*{GEX,VDJ}*R1.fq.gz; do
    pool_pattern=$(basename $file | grep -o P[1-9])
    lib_pattern=$(basename $file | grep -oE "GEX|VDJ")
    output_file="$OUTPUT/${pool_pattern}_${lib_pattern}.csv"
    almnt_file="$ALIGNMENT_DIR/${pool_pattern}_${lib_pattern}.out"
    jobname="Extract_Cellbarcode_${pool_pattern}_${lib_pattern}"
    echo "$file" "$output_file" "$almnt_file"
    sbatch --job-name="$jobname" \
           --error="${jobname}_err.txt" \
           --output="${jobname}_output.out" \
           --cpus-per-task=2 \
           --mem-per-cpu=4G \
           --wrap="python3 $SCRIPT -hmmer $almnt_file -R1 $file -out $output_file"
done