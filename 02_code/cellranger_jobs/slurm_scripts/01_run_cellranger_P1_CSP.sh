#!/bin/bash
#SBATCH --job-name=cellranger_P1_CSP
#SBATCH --error=car_t_sc/01_data/raw/cellranger_CSP_only/err.txt
#SBATCH --output=car_t_sc/01_data/raw/cellranger_CSP_only/output.out
#SBATCH --mail-type=ALL #if you want to receive emails when your job starts/fails/finishes, very handy
#SBATCH --mail-user=david.stumpf@stud-mail.uni-wuerzburg.de

##change working directory
cd /home/s377963/car_t_sc/01_data/raw/cellranger_CSP_only || exit

## Run the Cell Ranger command
/home/s377963/cellranger/cellranger-8.0.1/bin/cellranger count \
  --id CSP_p1_FB \
  --transcriptome /home/s377963/cellranger/refgenomes/refdata-gex-GRCm39-2024-A/ \
  --create-bam false \
  --feature-ref /home/s377963/car_t_sc/01_data/raw/cellranger/2024-06-07_24054SC_Luu_P1_cellranger/count/feature_reference.csv \
  --libraries /home/s377963/car_t_sc/02_code/cellranger_jobs/library_CSVs/CSP_only/P1.csv \
  --chemistry SC-FB
