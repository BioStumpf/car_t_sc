#!/bin/bash

#specify Data and output directories
DATA_DIR="$HOME/car_t_sc/01_data"
FASTQS_DIR="$DATA_DIR/raw/raw/2024-06-17_Maik_Luu_24054SC_raw_FASTQ/FASTQ_by_lib"
OUTPUT="$DATA_DIR/processed/count_cart_receptor"
R1_INPUT="$DATA_DIR/processed/bowtie2_mapped_carts"

for file in $FASTQS_DIR/*{GEX,VDJ}*R1.fq.gz; do
    pool_pattern=$(basename $file | grep -o P[1-9])
    lib_pattern=$(basename $file | grep -oE "GEX|VDJ")
    output_file="$OUTPUT/${pool_pattern}_${lib_pattern}.csv"
    R1_file="$R1_INPUT/${pool_pattern}_${lib_pattern}"
    # echo "$R1_file" "$output_file"
    sbatch --job_name=$
done
