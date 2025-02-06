#!/bin/bash

DATA_DIR="$HOME/car_t_sc/01_data"
INDEX="$DATA_DIR/processed/OVA_receptor/hmmer_hits_bowtie2/VDJ_Indx/vdj-GRCm38-alts-ensembl-7.0.0" #Genome_Indx/refdata-gex-GRCm39-2024-A or VDJ_Indx/vdj-GRCm38-alts-ensembl-7.0.0
FASTQS_DIR="$DATA_DIR/processed/OVA_receptor/hmmersearch_alignment_outputs/fastq_subsets"
OUTPUT="$DATA_DIR/processed/OVA_receptor/hmmer_hits_bowtie2/VDJ_bowtie2_output"

for file in "$FASTQS_DIR"/*.fq.gz; do #do VDJ or GEX
    pool_pattern=$(basename "$file" | grep -o P[1-9])
    lib_pattern=$(basename "$file" | grep -oE "GEX|VDJ")
    output_file="$OUTPUT/${pool_pattern}_${lib_pattern}.sam" 
    jobname="CART_mapping_${pool_pattern}_${lib_pattern}"
    # echo "$jobname" $pattern $output_file $file
    sbatch --job-name="$jobname" \
           --error="${jobname}_err.txt" \
           --output="${jobname}_output.out" \
           --cpus-per-task=10 \
           --mem-per-cpu=10G \
           --wrap="bowtie2 -x "$INDEX" -U "$file" -S "$output_file" --no-unal"
done


    # sbatch --job-name="Build_Indx" \
    #        --error="Build_Indx_err.txt" \
    #        --output="Build_Indx_output.out" \
    #        --cpus-per-task=25 \
    #        --mem-per-cpu=4G \
    #        --wrap="bowtie2-build cellranger/refgenomes/refdata-gex-GRCm39-2024-A/fasta/genome.fa car_t_sc/01_data/processed/OVA_receptor/hmmer_hits_bowtie2/Genome_Indx/refdata-gex-GRCm39-2024-A"

    #for vdj
    #  sbatch --job-name="bowtie2_build"        --error="bowtie2_build_err.txt"        --output="bowtie2_build_output.out"        --cpus-per-task=50        --mem=20G        --wrap="bowtie2-build $HOME/cellranger/refvdj/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0/fasta/regions.fa \
    #                     $HOME/car_t_sc/01_data/processed/OVA_receptor/hmmer_hits_bowtie2/VDJ_Indx/vdj-GRCm38-alts-ensembl-7.0.0"