#!/bin/bash

path_to_bam="$1"
output_file="$2"

samtools view $path_to_bam | awk -F'\t' '
                   {
                   raw = corrected = ""
                   for(i=1; i<=NF; i++) {
                   if ($i ~ /^CR:Z:/) { raw=substr($i,6) }
                   if ($i ~ /^CB:Z:/) { corrected=substr($i,6) }
                   }
                   if (raw && corrected) { print raw, corrected }
                   }' | sort | uniq > "$output_file"