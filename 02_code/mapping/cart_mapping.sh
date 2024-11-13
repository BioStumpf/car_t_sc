#!/bin/bash

DATA_DIR="$HOME/car_t_sc/01_data"
INDEX_DIR="$DATA_DIR/processed/mapping_index"
REFERENCE_DIR="$DATA_DIR/reference_cart"
FASTQS_DIR="$DATA_DIR/raw/raw/2024-06-17_Maik_Luu_24054SC_raw_FASTQ/FASTQ_by_lib"
OUTPUT="$DATA_DIR/processed/bowtie2_mapped_carts"
SLURM_OUTPUTS="$HOME/car_t_sc/02_code/mapping/ouputs_errors"

if ! ls "$INDEX_DIR/CD19_R11"*.bt2 1> /dev/null 2>&1; then
    echo "Index does not exist, creating Index first"
    bowtie2-build "$REFERENCE_DIR/CD19_R11.fa" "$INDEX_DIR/CD19_R11"
else 
    echo "Index exists"
fi

for file in "$FASTQS_DIR"/*VDJ*R2.fq.gz; do
    pattern=$(basename "$file" | grep -o P[1-9])
    output_file="$OUTPUT/${pattern}_VDJ.sam"
    jobname="CART_mapping_${pattern}"
    echo "Submitting ${jobname}"
    sbatch --job-name="$jobname" \
           --error="$SLURM_OUTPUTS/${jobname}_err.txt" \
           --output="$SLURM_OUTPUTS/${jobname}_output.out" \
           --cpus-per-task=2 \
            --mem-per-cpu=2G \
           --wrap="bowtie2 -x "$INDEX_DIR/CD19_R11" -U "$file" -S "$output_file" --no-unal"
            # cart_mapping_inner_func.sh "$INDEX_DIR/CD19_R11" "$file" "$output_file"
    echo "Successfully submitted ${jobname}"
done


# echo "$file"
# echo "$ouput_file"
# bowtie2 -x "$INDEX_DIR/CD19_R11" -U "$file" -S "$output_file" --no-unal
# bowtie2 -x "$INDEX_DIR/CD19_R11" -U "$FASTQS_DIR/24054SC_Luu_P1_D0_GEX_B11_TT_R2.fq.gz" -S "$OUTPUT/P1.sam" --no-unal
# bowtie2-build ~/car_t_sc/01_data/reference_cart/R11_scFV.fa ~/car_t_sc/01_data/processed/mapping_index/R11_scFV
# bowtie2 -x ~/car_t_sc/01_data/processed/mapping_index/R11_scFV -U  <(zcat ~/car_t_sc/01_data/raw/raw/2024-06-17_Maik_Luu_24054SC_raw_FASTQ/FASTQ_by_lib/24054SC_Luu_P1_D0_GEX_B11_TT_R2.fq.gz | head -n 1000) -S /home/s377963/car_t_sc/01_data/processed/bowtie2_mapped_carts/P1.sam --no-unal