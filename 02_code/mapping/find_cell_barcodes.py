#!/usr/bin/env python3

import argparse
import gzip
import numpy as np
import pandas as pd

class aligned_read_info:
    def __init__(self, id, motif, cellbarcode=None):
        self.id = id
        self.motif = motif
        self.cellbarcode = cellbarcode

def find_read_in_sam_reads(read_id, sam_aligned_reads):
    aligned_read_ids = [read.id for read in sam_aligned_reads]
    if read_id in aligned_read_ids:
        read_id_idx = aligned_read_ids.index(read_id)
    else:
        read_id_idx = -1
    return read_id_idx


def run(args):
    sam_file = open(args.sam)
    fastq_file = gzip.open(args.fastq, 'rt')
    fout = open(args.output, "w")

    #First we need the IDs from the SAM file, ergo we need to read it and store the IDs in a list
    aligned_reads = []
    for line in sam_file:
        #Skip the header since its not interesting here
        if line.startswith('@'):
            continue
        #split the line by tab and take the first entry == id
        line_elements = line.split('\t')
        id = line_elements[0]
        mapped_motif = line_elements[2]
        aligned_reads.append(aligned_read_info(id, mapped_motif))

    # fout.write('\n'.join(aligned_ids) + '\n')
    # fout.close()

    #second read the ids from the fastq file line by line and check if the id matches with the ids extracted from the sam file
    # aligned_id =''
    # read_indx = -1
    for line in fastq_file:
        if line.startswith('@'):
            #read only the first id part, since bowtie2 (the aligner that generated the SAM file) strips the id of the rest anyway
            read_id = line.split()[0][1:]
            #find the index of the current read_id within the sam file 
            read_indx = find_read_in_sam_reads(read_id, aligned_reads)
            #check if the current read_id present in the sam file
            # if read_indx != -1:
                #if yes, safe the id and the motif its binding to
                # aligned_id = read_id
                # motif = aligned_reads[read_indx].motif
            #after reading the id line, the sequence line is read, if you safed the id before, its because it matches the aligned ids from the SAM file, ergo we extract the 10x sequence
        elif line[0] != '@' and read_indx != -1:
            #extract the first 16 bases of the sequence (==10x seuquence)
            cell_barcode = line.strip()[:16]
            aligned_reads[read_indx].cellbarcode = cell_barcode
            #write the 10x barcode into the output file
            # fout.write(cell_barcode + '\t' + motif + '\n')
            #set the aligned_id to nothing to avoid reading none sequence lines (the + and the quality line)
            # aligned_id = ''
            read_indx = -1
            #after reading the sequence line, there are 2 more lines, these do not contain the '@' at the beginning of the line, however the aligned_id was previously set to nothing, so the elif condition from before is not met
        else:
            continue

    #need function that checks the aligned_reads list for multiplets, ergo counts how often specific 10x barcodes occur 
    #find unique motifs and cellbarcodes 
    unq_motifs = np.unique([read.motif for read in aligned_reads])
    unq_cellbarcodes = np.unique([read.cellbarcode for read in aligned_reads])
    #write all motifs into the output file
    motifs_print = ','.join(unq_motifs)
    fout.write(f'Cellbarcode,ReadIDs,{motifs_print}\n')
    #iterate through all barcodes
    for cellbarcode in unq_cellbarcodes:
        fout.write(cellbarcode)
        #subset the aligned reads to each unique barcode at a time
        subset = list(filter(lambda x: x.cellbarcode == cellbarcode, aligned_reads))
        #get the counts for each motif (whatever was mapped to)
        subset_motifs = pd.Index([read.motif for read in subset])
        IDs = ' '.join([read.id for read in subset])
        fout.write(f',{IDs}')
        motif_count = subset_motifs.value_counts()
        #iterate through each motif once 
        for motif in unq_motifs:
            #if the motif is within the motifs of the cellbarcode subset -> get its count and write it to the output file
            if motif in motif_count.index:
                motif_indx = motif_count.index.get_loc(motif)
                # motif_name = motif_count.index[motif_indx]
                motif_value = motif_count.values[motif_indx]
                # print(f'{motif_name}: {motif_value}')
            #if the motif is not within the motifs of the cellbarcode subset, set its cout to 0
            else: 
                motif_value = 0
            fout.write(f',{motif_value}')
        fout.write('\n')
    fout.close()

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