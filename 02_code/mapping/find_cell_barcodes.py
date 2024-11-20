#!/usr/bin/env python

import argparse
import gzip

def run(args):
    sam_file = open(args.sam)
    fastq_file = gzip.open(args.fastq, 'rt')
    fout = open(args.output, "w")

    #First we need the IDs from the SAM file, ergo we need to read it and store the IDs in a list
    aligned_ids = []
    aligned_motifs = []
    for line in sam_file:
        #Skip the header since its not interesting here
        if line.startswith('@'):
            continue
        #split the line by tab and take the first entry == id
        line_elements = line.split('\t')
        id = line_elements[0]
        mapped_motif = line_elements[2]
        aligned_ids.append(id)
        aligned_motifs.append(mapped_motif)
    # fout.write('\n'.join(aligned_ids) + '\n')
    # fout.close()
    #second read the ids from the fastq file line by line and check if the id matches with the ids extracted from the sam file
    aligned_id =''
    for line in fastq_file:
        if line.startswith('@'):
            #read only the first id part, since bowtie2 (the aligner that generated the SAM file) strips the id of the rest anyway
            read_id = line.split()[0][1:]
            # fout.write(read_id)
            # break
            #check if the current read_id present in the sam file
            if read_id in aligned_ids:
                #if yes, safe the id and the motif its binding to
                aligned_id = read_id
                position = aligned_ids.index(aligned_id)
                motif = aligned_motifs[position]
            #after reading the id line, the sequence line is read, if you safed the id before, its because it matches the aligned ids from the SAM file, ergo we extract the 10x sequence
        elif line[0] != '@' and aligned_id != '':
            #extract the first 16 bases of the sequence (==10x seuquence)
            cell_barcode = line.strip()[:16]
            #write the 10x barcode into the output file
            fout.write(cell_barcode + '\t' + motif + '\n')
            #set the aligned_id to nothing to avoid reading none sequence lines (the + and the quality line)
            aligned_id = ''
            #after reading the sequence line, there are 2 more lines, these do not contain the '@' at the beginning of the line, however the aligned_id was previously set to nothing, so the elif condition from before is not met
        else:
            continue
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