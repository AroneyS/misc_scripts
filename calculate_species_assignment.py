#!/usr/bin/env python3

"""
Author: Samuel Aroney
Summarise condense tables by calculating % species assignment
Sum of coverage values with species assignment divided by total coverage
"""

import sys
import os
import argparse
import logging
import polars as pl
from tqdm import tqdm

CONDENSED_COLUMNS = {
    "sample": str,
    "coverage": float,
    "taxonomy": str,
}

def pipeline(input_folder, output_file):
    logging.info(f"Collecting condensed tables from {input_folder}...")
    condensed_tables = []
    for root, _, files in os.walk(input_folder):
        for file in files:
            if file.endswith('condensed.csv'):
                condensed_tables.append(os.path.join(root, file))

    logging.info(f"Found {len(condensed_tables)} condensed tables")

    header = True
    with open(output_file, 'w') as f:
        for file in tqdm(condensed_tables):
            logging.debug(f"Processing {file}")
            table = (
                pl.scan_csv(file, separator="\t", dtypes=CONDENSED_COLUMNS)
                .with_columns(species = pl.col("taxonomy").str.contains("s__"))
                .group_by("sample", "species")
                .agg(pl.sum("coverage"))
                .collect(streaming=True)
                .pivot(index="sample", columns="species", values="coverage")
            )

            if "true" not in table.columns:
                table = table.with_columns(true = pl.lit(0))
            if "false" not in table.columns:
                table = table.with_columns(false = pl.lit(0))

            table = (
                table
                .with_columns(fraction = pl.col("true") / (pl.col("true") + pl.col("false")))
                .select("sample", "fraction")
            )

            if header:
                f.write("\t".join(table.columns) + "\n")
                header = False

            for row in table.iter_rows():
                f.write("\t".join(map(str, row)) + "\n")

    return


def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--debug", help="output debug information", action="store_true")
    parser.add_argument("--quiet", help="only output errors", action="store_true")

    parser.add_argument("--input", help="Input folder containing condense tables")
    parser.add_argument("--output", help="Output summary table")

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format="%(asctime)s %(levelname)s: %(message)s", datefmt="%Y/%m/%d %I:%M:%S %p")

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    pipeline(args.input, args.output)

    logging.info(f"Output written to {args.output}")
    logging.info("Done")

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
