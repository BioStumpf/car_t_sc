import argparse
import pandas as pd 

def run(args):
    with open(args.output, 'w') as out:
        with open(args.hmmer, 'r') as hmmer:
            # First we need the IDs from the hmmer output file
            for line in hmmer:
                # Skip the header since it's not interesting here
                if line.startswith('#'):
                    continue
                # Split the line by tab and take the first entry == id
                line_elements = line.split()
                read_id = line_elements[0]
                out.write(f'{read_id}\n')

def main():
    parser=argparse.ArgumentParser(description="Extract all IDs from hmmer mapping")
    parser.add_argument("-hmmer",help="hmmer search file" ,dest="hmmer", type=str, required=True)
    parser.add_argument("-out",help="txt output filename" ,dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()