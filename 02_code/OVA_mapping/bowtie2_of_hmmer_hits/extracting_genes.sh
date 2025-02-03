#!/bin/bash

GENES="$HOME/cellranger/refgenomes/refdata-gex-GRCm39-2024-A/genes/genes.gtf.gz"
INPT_DIR="$HOME/car_t_sc/01_data/processed/OVA_receptor/hmmer_hits_bowtie2/bowtie2_output"
FCOUNT_DIR="$HOME/car_t_sc/01_data/processed/OVA_receptor/hmmer_hits_bowtie2/feature_counts"
GENE_ID_DIR="$HOME/car_t_sc/01_data/processed/OVA_receptor/hmmer_hits_bowtie2/annotated_genes_for_bowtie_hits"


for file in $INPT_DIR/*.sam; do
    pool_pattern=$(basename $file | grep -o "P[1-9]")
    lib_pattern=$(basename $file | grep -oE "GEX|VDJ")
    fcounts="${FCOUNT_DIR}/${pool_pattern}_${lib_pattern}.txt"
    fcounts_bed="${FCOUNT_DIR}/${pool_pattern}_${lib_pattern}.bed"
    gene_id_annotation="${GENE_ID_DIR}/${pool_pattern}_${lib_pattern}.bed"
    zcat $GENES | featureCounts -a /dev/stdin -o $fcounts -t gene -g gene_name $file
    awk '$NF > 0' "$fcounts" > "${fcounts}.filtered"
    awk 'BEGIN {OFS="\t"} {if ($1 != "GeneID") print $2, $3-1, $4, $1, ".", $5}' "${fcounts}.filtered" > $fcounts_bed
    samtools view -Sb $file | bedtools intersect -a /dev/stdin -b $fcounts_bed -wa -wb -bed> $gene_id_annotation
    rm -r $FCOUNT_DIR

done


    # awk '$NF > 0' "$output_file" > "${output_file}.filtered"

# zcat /home/s377963/cellranger/refgenomes/refdata-gex-GRCm39-2024-A/genes/genes.gtf.gz | \
# featureCounts -a /dev/stdin -o /home/s377963/car_t_sc/01_data/processed/OVA_receptor/hmmer_hits_bowtie2/feature_counts/P1_GEX_gene_counts.txt \
# -t gene -g gene_name /home/s377963/car_t_sc/01_data/processed/OVA_receptor/hmmer_hits_bowtie2/bowtie2_output/P1_GEX.sam