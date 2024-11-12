#!/bin/bash
DIR="$HOME/car_t_sc/01_data/processed/mapping_index"

if ! ls "$DIR/CD19_R11"*.bt2 1> /dev/null 2>&1; then
    echo "Index does not exist, creating Index first"
    bowtie2-build "$HOME/car_t_sc/01_data/reference_cart/CD19_R11.fa" "$HOME/car_t_sc/01_data/processed/mapping_index/CD19_R11"
else 
    echo "Index exists"
fi

#~/car_t_sc/01_data/processed/mapping_index/R11_scFV
# bowtie2-build ~/car_t_sc/01_data/reference_cart/CD19_R11.fa ~/car_t_sc/01_data/processed/mapping_index/R11_scFV
# bowtie2 -x ~/car_t_sc/01_data/processed/mapping_index/R11_scFV -U ~/car_t_sc/01_data/raw/raw/2024-06-17_Maik_Luu_24054SC_raw_FASTQ/FASTQ_by_lib/24054SC_Luu_P1_D0_GEX_B11_TT_R2.fq.gz -S /home/s377963/car_t_sc/01_data/processed/bowtie2_mapped_carts/P1.sam --no-unal


# bowtie2-build ~/car_t_sc/01_data/reference_cart/R11_scFV.fa ~/car_t_sc/01_data/processed/mapping_index/R11_scFV
# bowtie2 -x ~/car_t_sc/01_data/processed/mapping_index/R11_scFV -U  <(zcat ~/car_t_sc/01_data/raw/raw/2024-06-17_Maik_Luu_24054SC_raw_FASTQ/FASTQ_by_lib/24054SC_Luu_P1_D0_GEX_B11_TT_R2.fq.gz | head -n 1000) -S /home/s377963/car_t_sc/01_data/processed/bowtie2_mapped_carts/P1.sam --no-unal