#!/bin/bash

cellranger_outs_dir="$HOME/car_t_sc/01_data/raw/cellranger_multi_CAR" 
correction_folder="$HOME/car_t_sc/01_data/processed/cellbarcode_correction"

for folder in $cellranger_outs_dir/*; do
    pattern=$(basename $folder) 
    path_to_bam="$cellranger_outs_dir/${pattern}/${pattern}/outs/per_sample_outs/${pattern}/count/sample_alignments.bam"
    jobname="${pattern}_brcd_correction"
    output_file="${correction_folder}/${pattern}_barcode_mapping.tsv"
    # echo $jobname $pattern $path_to_bam
    sbatch --job-name=$jobname \
           --error="${jobname}_err.txt" \
           --output="${jobname}.out" \
           --mem=20G \
           ./process_bam.sh $path_to_bam $output_file
done

# --wrap="samtools view $path_to_bam | awk -F'\t' '{
#                    raw = corrected = ""
#                    for(i=1; i<=NF; i++) {
#                    if ($i ~ /^CR:Z:/) { raw=substr($i,6) }
#                    if ($i ~ /^CB:Z:/) { corrected=substr($i,6) }
#                    }
#                    if (raw && corrected) { print raw, corrected }
#                    }' | sort | uniq > "${correction_folder}/${pattern}_barcode_mapping.tsv""

   
    # if [[ $pattern == "P4" || $pattern == "P7"]]; then 
    #     continue
    # fi


    # samtools view $cellranger_outs_dir | awk -F'\t' '{
    #     raw = corrected = ""
    #     for(i=1; i<=NF; i++) {
    #         if ($i ~ /^CR:Z:/) { raw=substr($i,6) }
    #         if ($i ~ /^CB:Z:/) { corrected=substr($i,6) }
    #     }
    #     if (raw && corrected) { print raw, corrected }
    # }' | sort | uniq > "${correction_folder}/${pattern}_barcode_correction.tsv"