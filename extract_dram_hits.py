#!/usr/bin/env python3

"""
Author: Samuel Aroney
Extract DRAM-matched sequences from fasta file by name and coordinates

Test:
python extract_dram_hits.py --input test/test_contigs.fa --dram test/test_dram.tsv --output test/test_extract_output.fna
python extract_dram_hits.py --input test/test_contigs.fa --dram test/test_dram.tsv --output test/test_extract_output.faa --translate
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
    parser.add_argument("--translate", help="Translate sequences to protein", action="store_true")

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    logging.info("Reading input fasta file")
    sequences = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))

    def extract_sequence(contig, begin, end, strand, sequences=sequences, translate=args.translate):
        sequence = sequences[contig].seq[begin - 1:end]
        if strand == -1:
            sequence = sequence.reverse_complement()

        if translate:
            sequence = sequence.translate()

        return str(sequence)

    logging.info("Processing DRAM file")
    output = (
        pl.read_csv(args.dram, separator="\t")
        .rename({"": "name"})
        .with_columns(
            contig = pl.concat_str(["fasta", "scaffold"], separator="_"),
        )
        .select(
            "name",
            info = pl.struct(["contig", "start_position", "end_position", "strandedness"])
                .apply(lambda x: extract_sequence(x["contig"], x["start_position"], x["end_position"], x["strandedness"]))
            )
    )

    logging.info("Writing output fasta file")
    with open(args.output, "w") as f:
        for name, seq in output.iter_rows():
            f.write(">%s\n" % name)
            f.write("%s\n" % seq)

    logging.info("Done")

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
