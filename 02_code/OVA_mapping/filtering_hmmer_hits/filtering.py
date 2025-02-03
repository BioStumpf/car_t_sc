import argparse
import pandas as pd 
import re

def run(args):
    to_filter = []
    with open(args.annotated, 'r') as antd:
        for line in antd:
            line_elements = line.split()
            tcr_name = re.compile(r"^Tr")
            if tcr_name.match(line_elements[15]):
                continue
            to_filter.append(line_elements[3])

    with open(args.output, 'w') as out:
        with open(args.hmmer, 'r') as hmmer:
            for line in hmmer:
                if line.startswith('#'):
                    continue
                line_elements = line.split()
                id = line_elements[0]
                if id not in to_filter:
                    out.write(line)

def main():
    parser=argparse.ArgumentParser(description="Extract all IDs from hmmer mapping")
    parser.add_argument("-hmmer_annotated",help="hmmer hits annotated to the genome" ,dest="annotated", type=str, required=True)
    parser.add_argument("-hmmer",help="hmmer search file" ,dest="hmmer", type=str, required=True)
    parser.add_argument("-out",help="txt output filename" ,dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()