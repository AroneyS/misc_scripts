#!/usr/bin/env python3

"""
Author: Samuel Aroney
Choose quality bins from CheckM2 output
"""

import sys
import os
import argparse
import logging
import polars as pl

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument('--checkm2-output', help='CheckM2 quality_report.tsv', required=True)
    parser.add_argument("--bins", help="Bins list", required=True)
    parser.add_argument("--output", help="Output file", required=True)
    parser.add_argument("--min-completeness", help="Minimum bin completeness", type=int, default=70)
    parser.add_argument("--max-contamination", help="Maximum bin contamination", type=int, default=10)

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

    bins = (
        pl.read_csv(args.bins, has_header=False, new_columns=["path"])
        .with_columns(
            Name = pl.col("path").apply(lambda x: os.path.splitext(os.path.basename(x))[0])
            )
    )

    output = (
        pl.read_csv(args.checkm2_output, separator="\t")
        .join(bins, on="Name")
        .filter(pl.col("Completeness") >= args.min_completeness)
        .filter(pl.col("Contamination") <= args.max_contamination)
        .select("path")
    )

    output.write_csv(args.output, separator = "\t", has_header=False)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
