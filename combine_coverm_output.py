#!/usr/bin/env python3

"""
Author: Samuel Aroney
Combine output files from CoverM across samples
"""

import sys
import argparse
import logging
import polars as pl

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument("--coverm-list", help="File containing list of CoverM outputs to combine", required=True)
    parser.add_argument("--output", help="Output combined file", required=True)

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

    with open(args.coverm_list) as f:
        coverm_files = [line.strip() for line in f]

    first = True
    for file in coverm_files:
        if first:
            first = False
            df = pl.scan_csv(file, separator="\t")
        else:
            df = df.join(pl.scan_csv(file, separator="\t"), on="Genome", how="inner")

    df.collect().write_csv(args.output, separator="\t")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
