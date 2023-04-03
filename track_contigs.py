#!/usr/bin/env python3

"""
Author: Samuel Aroney
Track contigs from assembly to genomes by name
"""

import os
import sys
import argparse
import logging
from Bio import SeqIO
import polars as pl


def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument("--assembly", help="Assembly file")
    parser.add_argument("--genomes", nargs='+', help="Genome file/s")
    parser.add_argument("-o", "--outfile", help="Output file")

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

    contigs = []
    with open(args.assembly) as f:
        for record in SeqIO.parse(f, "fasta"):
            contigs.append(record.id)
    assembly_contigs = pl.DataFrame({
        "assemblies": os.path.splitext(os.path.basename(args.assembly))[0],
        "contigs": contigs,
    })

    genome_contigs = pl.DataFrame({
        "genome_paths": args.genomes,
    }).select(
        pl.col("genome_paths").apply(lambda x: [r.id for r in SeqIO.parse(x, "fasta")]).alias("contigs"),
        pl.col("genome_paths").apply(lambda x: os.path.splitext(os.path.basename(x))[0]).alias("genomes"),
    ).explode("contigs")

    matches = assembly_contigs.join(genome_contigs, on="contigs", how="left")
    matches.write_csv(args.outfile, sep="\t")

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
