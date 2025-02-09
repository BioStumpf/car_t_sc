#!/bin/bash

csv_dir="$HOME/car_t_sc/02_code/cellranger_jobs/multi_config_csvs"

for config in $csv_dir/*; do
  pool_pattern=$(basename $config | grep -o p[1-9]| tr '[:lower:]' '[:upper:]')
  # Skip P1 (already processed)
  if [[ $pool_pattern == "P4" || $pool_pattern == "P7" ]]; then
    jobname="${pool_pattern}_multi"
    echo $HOME/car_t_sc/01_data/raw/cellranger_multi_CAR/${pool_pattern}
    sbatch --job-name=$jobname \
           --error="${jobname}_err.txt" \
           --output="${jobname}.out" \
           --cpus-per-task=16 \
           --mem=100G \
           --wrap="cd $HOME/car_t_sc/01_data/raw/cellranger_multi_CAR/${pool_pattern} && cellranger multi --id $pool_pattern --csv $config"
  fi
done

