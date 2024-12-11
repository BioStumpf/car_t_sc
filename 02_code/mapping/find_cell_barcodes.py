#!/usr/bin/env python3

import argparse
import gzip
import numpy as np
import pandas as pd

# class aligned_read_info:
#     def __init__(self, id, motif, cellbarcode):
#         self.id = id
#         self.motif = motif
#         self.cellbarcode = cellbarcode

def run(args):
    sam_file = open(args.sam)
    fastq_file = gzip.open(args.fastq, 'rt')
    #prepare fastqs for finding the matched read id and corresponding cellbarcodes
    fastq_lines = fastq_file.readlines()
    stripped_ids_dict = {line.split()[0]: idx for idx, line in enumerate(fastq_lines[0:len(fastq_lines):4])} #use dictionary for o(1) searching
    #First we need the IDs from the SAM file, ergo we need to read it and store the IDs in a list
    rows = []
    for line in sam_file:
        #Skip the header since its not interesting here
        if line.startswith('@'):
            continue
        #split the line by tab and take the first entry == id
        line_elements = line.split('\t')
        read_id = f'@{line_elements[0]}'
        mapped_motif = line_elements[2]
        #now search the mapped read_id within the fastq file
        read_indx = stripped_ids_dict[read_id]
        #look for the barcode index (note a fastq file is a multiple of 4)
        barcode_indx = (read_indx * 4 ) + 1
        barcode = fastq_lines[barcode_indx][:16]
        rows.append([barcode, read_id, mapped_motif])
    sam_file.close()
    fastq_file.close()

    aligned_reads = pd.DataFrame(rows, columns=["Cellbarcode", "ReadIDs", "Motif"])
    reads_grouped = aligned_reads.groupby("Cellbarcode").agg({"ReadIDs": lambda x: "\t".join(map(str, x))})
    reads_pivot = aligned_reads.pivot_table(index="Cellbarcode", columns="Motif", aggfunc="size", fill_value=0)
    motif_count_table = pd.merge(reads_grouped, reads_pivot, on="Cellbarcode")
    motif_count_table.to_csv(args.output)

        # aligned_reads.append(aligned_read_info(read_id, mapped_motif, barcode))
        # aligned_reads.loc[len(aligned_reads)] = [barcode, read_id, mapped_motif]
    #create a pandas dataframe to append the information to
    # aligned_reads = pd.DataFrame(columns=["Cellbarcode", "ReadIDs", "Motif"])
    # fout = open(args.output, "w")
    


    #need function that checks the aligned_reads list for multiplets, ergo counts how often specific 10x barcodes occur 
    #find unique motifs and cellbarcodes 
    # unq_motifs = np.unique([read.motif for read in aligned_reads])
    # unq_cellbarcodes = np.unique([read.cellbarcode for read in aligned_reads])
    # #write all motifs into the output file
    # motifs_print = ','.join(unq_motifs)
    # fout.write(f'Cellbarcode,ReadIDs,{motifs_print}\n')
    # #iterate through all barcodes
    # for cellbarcode in unq_cellbarcodes:
    #     fout.write(cellbarcode)
    #     #subset the aligned reads to each unique barcode at a time
    #     subset = list(filter(lambda x: x.cellbarcode == cellbarcode, aligned_reads))
    #     #get the counts for each motif (whatever was mapped to)
    #     subset_motifs = pd.Index([read.motif for read in subset])
    #     IDs = '\t'.join([read.id for read in subset])
    #     fout.write(f',{IDs}')
    #     motif_count = subset_motifs.value_counts()
    #     #iterate through each motif once 
    #     for motif in unq_motifs:
    #         #if the motif is within the motifs of the cellbarcode subset -> get its count and write it to the output file
    #         if motif in motif_count.index:
    #             motif_indx = motif_count.index.get_loc(motif)
    #             # motif_name = motif_count.index[motif_indx]
    #             motif_value = motif_count.values[motif_indx]
    #             # print(f'{motif_name}: {motif_value}')
    #         #if the motif is not within the motifs of the cellbarcode subset, set its cout to 0
    #         else: 
    #             motif_value = 0
    #         fout.write(f',{motif_value}')
    #     fout.write('\n')
    # fout.close()

def main():
    parser=argparse.ArgumentParser(description="Find Fastq identity from SAM file in R1 of 10x genomics output. Output is the 10x sequence, ergo the cellular identity.")
    parser.add_argument("-sam",help="sequence alignment map file" ,dest="sam", type=str, required=True)
    parser.add_argument("-R1",help="R1 fastq input file" ,dest="fastq", type=str, required=True)
    parser.add_argument("-out",help="fastq output filename" ,dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()