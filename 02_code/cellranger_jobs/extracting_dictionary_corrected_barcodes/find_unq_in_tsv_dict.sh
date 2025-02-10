#!/bin/bash

non_unq_tsv="$HOME/car_t_sc/01_data/processed/cellbarcode_correction/raw_non_unq_barcodes"
output_dir="$HOME/car_t_sc/01_data/processed/cellbarcode_correction/unq_barcodes"

for file in $non_unq_tsv/*; do
    output_file="${output_dir}/$(basename $file)"
    sbatch --job-name=$jobname \
       --error="${jobname}_err.txt" \
       --output="${jobname}.out" \
       --mem=20G \
       --wrap="tsv-uniq -f 1 $file > $output_file"
done