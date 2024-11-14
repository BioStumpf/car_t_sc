#!/bin/bash

#check if we  have 4 arguments given, if not quit the job immediately
if [ $# -ne 3 ]; then 
    echo "Usage: $0 path_to_idx file_to_map output_file"
    exit 1
fi

#assign the variables for bowtie2 given by the input
idx_dir=$1
file=$2
output_file=$3

#run bowtie2
bowtie2 -x "$idx_dir/CD19_R11" -U "$file" -S "$output_file" --no-unal