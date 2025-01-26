#!/usr/bin/env python3

import argparse
import gzip
import pandas as pd 

def run(args):
    hmmer = open(args.hmmer)
    fastq_file = gzip.open(args.fastq, 'rt')
    #prepare fastqs for finding the matched read id and corresponding cellbarcodes
    fastq_lines = fastq_file.readlines()
    stripped_ids_dict = {line.split()[0]: idx for idx, line in enumerate(fastq_lines[0:len(fastq_lines):4])} #use dictionary for o(1) searching
    #First we need the IDs from the hmmer output file, ergo we need to read it and store the IDs in a list
    rows = []
    for line in hmmer:
        #Skip the header since its not interesting here
        if line.startswith('#'):
            continue
        #split the line by tab and take the first entry == id
        line_elements = line.split()
        read_id = f'@{line_elements[0]}'
        mapped_motif = line_elements[2]
        #now search the mapped read_id within the fastq file
        read_indx = stripped_ids_dict[read_id]
        #look for the barcode index (note a fastq file is a multiple of 4)
        barcode_indx = (read_indx * 4 ) + 1
        barcode = fastq_lines[barcode_indx][:16]
        rows.append([barcode, read_id, mapped_motif])
    hmmer.close()
    fastq_file.close()

    aligned_reads = pd.DataFrame(rows, columns=["Cellbarcode", "ReadIDs", "Motif"])
    reads_grouped = aligned_reads.groupby("Cellbarcode").agg({"ReadIDs": lambda x: "\t".join(map(str, x))})
    reads_pivot = aligned_reads.pivot_table(index="Cellbarcode", columns="Motif", aggfunc="size", fill_value=0)
    motif_count_table = pd.merge(reads_grouped, reads_pivot, on="Cellbarcode")
    motif_count_table.to_csv(args.output)

def main():
    parser=argparse.ArgumentParser(description="Find Fastq identity from SAM file in R1 of 10x genomics output. Output is the 10x sequence, ergo the cellular identity.")
    parser.add_argument("-hmmer",help="hmmer search file" ,dest="hmmer", type=str, required=True)
    parser.add_argument("-R1",help="R1 fastq input file" ,dest="fastq", type=str, required=True)
    parser.add_argument("-out",help="fastq output filename" ,dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()