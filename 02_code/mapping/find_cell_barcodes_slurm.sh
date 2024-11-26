#!/bin/bash

#specify Data and output directories
DATA_DIR="$HOME/car_t_sc/01_data"
FASTQS_DIR="$DATA_DIR/raw/raw/2024-06-17_Maik_Luu_24054SC_raw_FASTQ/FASTQ_by_lib"
OUTPUT="$DATA_DIR/processed/count_cart_receptor"
SAM_INPUT="$DATA_DIR/processed/bowtie2_mapped_carts"
SCRIPT_DIR="$HOME/car_t_sc/02_code/mapping"

for file in $FASTQS_DIR/*{GEX,VDJ}*R1.fq.gz; do
    pool_pattern=$(basename $file | grep -o P[1-9])
    lib_pattern=$(basename $file | grep -oE "GEX|VDJ")
    output_file="$OUTPUT/${pool_pattern}_${lib_pattern}.csv"
    sam_file="$SAM_INPUT/${pool_pattern}_${lib_pattern}.sam"
    jobname="Extract_Cellbarcode_${pool_pattern}_${lib_pattern}"
    echo "$file" "$output_file" "$sam_file"
    sbatch --job-name="$jobname" \
           --error="${jobname}_err.txt" \
           --output="${jobname}_output.out" \
           --cpus-per-task=2 \
           --mem-per-cpu=2G \
           --wrap="python3 ${SCRIPT_DIR}/find_cell_barcodes.py -sam '$sam_file' -R1 "$file" -out '$output_file'"
done
