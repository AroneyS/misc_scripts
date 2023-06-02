#!/usr/bin/env python3

"""
Author: Samuel Aroney
Extract DRAM-matched sequences from fasta file by name and coordinates

Test: python extract_dram_hits.py --input test/test_contigs.fa --dram test/test_dram.tsv --output test/test_extract_output.fa
"""

import sys
import argparse
import logging
from Bio import SeqIO
import polars as pl


def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument("--input", help="Input fasta file", required=True)
    parser.add_argument("--dram", help="Input DRAM tsv file", required=True)
    parser.add_argument("--output", help="Output fasta file", required=True)

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

    def extract_sequence(contig, begin, end, strand, sequences=sequences):
        sequence = sequences[contig].seq[begin - 1:end]
        if strand == -1:
            sequence = sequence.reverse_complement()

        return str(sequence)

    output = (
        pl.read_csv(args.dram, separator="\t")
        .rename({"": "name"})
        .select(
            "name",
            info = pl.struct(["scaffold", "start_position", "end_position", "strandedness"])
                .apply(lambda x: extract_sequence(x["scaffold"], x["start_position"], x["end_position"], x["strandedness"]))
            )
    )

    with open(args.output, "w") as f:
        for name, seq in output.iter_rows():
            f.write(">%s\n" % name)
            f.write("%s\n" % seq)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
