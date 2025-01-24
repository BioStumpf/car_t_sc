#!/usr/bin/env python3

from Bio.Data import CodonTable
from itertools import product
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
import gzip
import argparse

#helper function for avoiding errrors due to nucleotide sequences not being a multiple of 3
def multiple3(self): 
    rem = len(self) % 3
    if rem != 0:
        self = self[:-rem]
    return self
MutableSeq.multiple3 = multiple3

#main funtion capable of converting a nucleotide fastq.gz file into a fasta output, note each read_id is transferred 3 times, once for each open reading frame
#also, nucleotide triplets not encoding for any aminoacid at all are transferred into an arbitrary nucleotide, denoted by 'X'
def translate_fastq_to_fasta_out(args):
    output = args.output
    input = args.input
    with gzip.open(output, 'wt') as output_handle: 
        with gzip.open(input, "rt") as input_file:
            for record in SeqIO.parse(input_file, 'fastq'):
                seq = MutableSeq(record.seq)
                for start in range(0,3):
                    translated = seq[start:].multiple3().translate().replace('*', 'X')
                    new_record = SeqRecord(translated, id = record.id, description=f"frame{start+1}")
                    SeqIO.write(new_record, output_handle, "fasta")

def main():
    parser=argparse.ArgumentParser(description="Convert the fq.gz illumina sequencing fastqs which are nucleotide sequences into protein sequences, base on all 3 ORFs")
    parser.add_argument("-in",help=".fq.gz file containing the reads" ,dest="input", type=str, required=True)
    parser.add_argument("-out",help=".fa output filename which will contain protein sequences" ,dest="output", type=str, required=True)
    parser.set_defaults(func=translate_fastq_to_fasta_out)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()                    