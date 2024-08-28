#!/bin/bash
#SBATCH --job-name=cellranger_P1_multi
#SBATCH --error=car_t_sc/01_data/raw/cellranger_CSP_GEX/err.txt
#SBATCH --output=car_t_sc/01_data/raw/cellranger_CSP_GEX/output.out

##change working directory
cd /home/s377963/car_t_sc/01_data/raw/multi_test || exit

## Run the Cell Ranger command
/home/s377963/cellranger/cellranger-8.0.1/bin/cellranger multi \
  --id p1 \
  --csv /home/s377963/car_t_sc/02_code/cellranger_jobs/multi_config_csvs/p1_multi_config.csv
