#!/bin/bash

DATA_DIR="$HOME/car_t_sc/01_data"
MODEL_DIR="$DATA_DIR/processed/OVA_receptor/hmmer_models"
FASTQS_DIR="$DATA_DIR/processed/OVA_receptor/raw_FASTQS_protein"
OUTPUT_TBL="$DATA_DIR/processed/OVA_receptor/hmmersearch_alignment_outputs"
OUTPUT_TXT="$DATA_DIR/processed/OVA_receptor/hmmersearch_alignment_outputs/summary_outputs"
MODEL="$MODEL_DIR/TCRs.hmm"

for file in $FASTQS_DIR/*{GEX,VDJ}*.fa; do #do VDJ or GEX
    pool_pattern=$(basename "$file" | grep -o P[1-9]) 
    lib_pattern=$(basename "$file" | grep -oE "GEX|VDJ")
    output_txt="$OUTPUT_TXT/${pool_pattern}_${lib_pattern}.txt" #do VDJ or GEX
    output_table="$OUTPUT_TBL/${pool_pattern}_${lib_pattern}.out"
    jobname="OVA_mapping_${pool_pattern}_${lib_pattern}"
    echo "$jobname" "$file" "$pool_pattern" "$lib_pattern"
    sbatch --job-name="$jobname" \
           --error="${jobname}_err.txt" \
           --output="${jobname}_output.out" \
           --cpus-per-task=2 \
           --mem-per-cpu=2G \
           --wrap="hmmsearch -o $output_txt --tblout $output_table $MODEL $file"
done
