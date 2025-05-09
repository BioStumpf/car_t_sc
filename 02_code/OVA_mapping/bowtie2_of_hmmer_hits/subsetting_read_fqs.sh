#!/bin/bash

#specify Data and output directories
DATA_DIR="$HOME/car_t_sc/01_data"
FASTQS_DIR="$DATA_DIR/raw/raw/2024-06-17_Maik_Luu_24054SC_raw_FASTQ/FASTQ_by_lib"
OUTPUT="$DATA_DIR/processed/OVA_receptor/hmmersearch_alignment_outputs/fastq_subsets/R1"
HMMER_DIR="$DATA_DIR/processed/OVA_receptor/hmmersearch_alignment_outputs/ids_txts_only"

for file in $FASTQS_DIR/*{GEX,VDJ}*R1.fq.gz; do
    pool_pattern=$(basename $file | grep -o P[1-9])
    lib_pattern=$(basename $file | grep -oE "GEX|VDJ")
    output_file="$OUTPUT/${pool_pattern}_${lib_pattern}_R1.fq.gz"
    hmmer_hit="$HMMER_DIR/${pool_pattern}_${lib_pattern}.txt"
    jobname="Subset_fqs_${pool_pattern}_${lib_pattern}"
    echo "$file" "$output_file"
    sbatch --job-name="$jobname" \
           --error="${jobname}_err.txt" \
           --output="${jobname}_output.out" \
           --cpus-per-task=10\
           --mem-per-cpu=10G \
           --wrap="seqkit grep -f $hmmer_hit $file -o $output_file"
done