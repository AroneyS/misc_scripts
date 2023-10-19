#!/usr/bin/env python3

"""
Author: Samuel Aroney
Find lengths for each entry in fasta file
"""

import sys
import os
import argparse
from Bio import SeqIO

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input-file', help='List of fasta files', required=True)

    args = parser.parse_args(arguments)

    with open(args.input_file) as input_file:
        files = input_file.read().splitlines()

    print("\t".join(["fasta", "record", "contig_length"]))
    for file in files:
        filename = os.path.splitext(os.path.basename(file))[0]
        with open(file) as f:
            for record in SeqIO.parse(f, 'fasta'):
                print("\t".join([filename, record.id, str(len(record.seq))]))

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
