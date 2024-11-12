#!/bin/bash
#SBATCH --job-name=mapping_pools_to_CART
#SBATCH --error=err.txt
#SBATCH --output=output.out

DATA_DIR="$HOME/car_t_sc/01_data"
INDEX_DIR="$DATA_DIR/processed/mapping_index"
REFERENCE_DIR="$DATA_DIR/reference_cart"
FASTQS_DIR="$DATA_DIR/raw/raw/2024-06-17_Maik_Luu_24054SC_raw_FASTQ/FASTQ_by_lib"
OUTPUT="$DATA_DIR/processed/bowtie2_mapped_carts"

if ! ls "$INDEX_DIR/CD19_R11"*.bt2 1> /dev/null 2>&1; then
    echo "Index does not exist, creating Index first"
    bowtie2-build "$REFERENCE_DIR/CD19_R11.fa" "$INDEX_DIR/CD19_R11"
else 
    echo "Index exists"
fi

for file in "$FASTQS_DIR"/*GEX*R2.fq.gz; do
    pattern=$(basename "$file" | grep -o P[0-9])
    output_file="$OUTPUT/${pattern}_GEX.sam"
    bowtie2 -x "$INDEX_DIR/CD19_R11" -U "$file" -S "$output_file" --no-unal
done


# echo "$file"
# echo "$ouput_file"
# bowtie2 -x "$INDEX_DIR/CD19_R11" -U "$file" -S "$output_file" --no-unal
# bowtie2 -x "$INDEX_DIR/CD19_R11" -U "$FASTQS_DIR/24054SC_Luu_P1_D0_GEX_B11_TT_R2.fq.gz" -S "$OUTPUT/P1.sam" --no-unal
# bowtie2-build ~/car_t_sc/01_data/reference_cart/R11_scFV.fa ~/car_t_sc/01_data/processed/mapping_index/R11_scFV
# bowtie2 -x ~/car_t_sc/01_data/processed/mapping_index/R11_scFV -U  <(zcat ~/car_t_sc/01_data/raw/raw/2024-06-17_Maik_Luu_24054SC_raw_FASTQ/FASTQ_by_lib/24054SC_Luu_P1_D0_GEX_B11_TT_R2.fq.gz | head -n 1000) -S /home/s377963/car_t_sc/01_data/processed/bowtie2_mapped_carts/P1.sam --no-unal