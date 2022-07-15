#!/usr/bin/env python3

"""
Author: Samuel Aroney
Extract rRNA sequences from fasta file by name and coordinates
"""

import os
import sys
import argparse
import logging
from Bio import SeqIO
import pandas as pd


def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument("--input", type=argparse.FileType("r"), help="Input fasta file")
    parser.add_argument("--dram-rrna", type=argparse.FileType("r"), help="Input DRAM rRNA tsv file")
    parser.add_argument("-o", "--outfile", type=argparse.FileType("w"), default=sys.stdout, help="Output fasta file")

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    sequences = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    rRNAs = pd.read_csv(args.dram_rrna, sep="\t")

    def extract_rRNAs(row, sequences):
        if row["strand"] == "-":
            sequence = sequences[row["scaffold"]].seq[row["begin"] - 1:row["end"]].reverse_complement()
        else:
            sequence = sequences[row["scaffold"]].seq[row["begin"] - 1:row["end"]]
        return str(sequence)

    rRNAs["sequence"] = rRNAs.apply(extract_rRNAs, sequences=sequences, axis=1)
    rRNAs.to_csv(args.outfile, sep="\t", index=False)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
