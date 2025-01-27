#!/bin/bash

#specify Data and output directories
FASTQS_DIR="$HOME/car_t_sc/01_data/processed/OVA_receptor/raw_FASTQS_protein"

for file in $FASTQS_DIR/*; do
    filename=$(basename $file)
    echo "$file"
    sbatch --job-name="unzip_$filename" \
           --error="unzip_${filename}_err.txt" \
           --output="unzip_${filename}_output.out" \
           --cpus-per-task=2 \
           --mem-per-cpu=2G \
           --wrap="gunzip $file"
done