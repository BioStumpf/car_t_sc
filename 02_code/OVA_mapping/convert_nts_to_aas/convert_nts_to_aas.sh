#!/bin/bash

#specify Data and output directories
DATA_DIR="$HOME/car_t_sc/01_data"
FASTQS_DIR="$DATA_DIR/raw/raw/2024-06-17_Maik_Luu_24054SC_raw_FASTQ/FASTQ_by_lib"
OUTPUT="$DATA_DIR/processed/OVA_receptor/raw_FASTQS_protein"
SCRIPT="$HOME/car_t_sc/02_code/OVA_mapping/convert_nts_to_aas/convert_nts_to_aas.py"


for file in $FASTQS_DIR/*{GEX,VDJ}*R2.fq.gz; do
    pool_pattern=$(basename $file | grep -o P[1-9])
    lib_pattern=$(basename $file | grep -oE "GEX|VDJ")
    output_file="$OUTPUT/${pool_pattern}_${lib_pattern}.fa"
    jobname="Extract_Cellbarcode_${pool_pattern}_${lib_pattern}"
    echo "$file" "$output_file"
    sbatch --job-name="$jobname" \
           --error="${jobname}_err.txt" \
           --output="${jobname}_output.out" \
           --cpus-per-task=10\
           --mem-per-cpu=2G \
           --wrap="seqkit translate -f 1,2,3 --clean '$file' -o '$output_file'"
done

        #    --wrap="python3 ${SCRIPT} -in '$file' -out '$output_file' -cpus 35"
