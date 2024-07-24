#!/usr/bin/env python3

"""
Author: Samuel Aroney
Filter Ibis proposed coassemblies to top N targets
"""

import os
import sys
import argparse
import logging
import polars as pl
import shutil


def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument("--naming", help="naming tsv with name and SampleID__", required=True)
    parser.add_argument("--jgi", help="JGI tsv with label (name) and filename", required=True)
    parser.add_argument("--input-dir", help="Directory containing downloads", required=True)
    parser.add_argument("--output-dir", help="Directory to write renamed files to", required=True)

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    naming = pl.read_csv(args.naming, separator="\t")
    jgi = pl.read_csv(args.jgi, separator="\t")
    conversion = (
        jgi
        .join(naming, left_on="label", right_on="name")
        .select("filename", "SampleID__")
    )
    logging.info(f"Processing {jgi.height} JGI samples into {conversion.height} renamed samples")

    samples = (
        pl.DataFrame({
            "filepath": [os.path.join(dp, f) for dp, dn, filenames in os.walk(args.input_dir) for f in filenames if f in conversion["filename"].to_list()]
            })
        .with_columns(filename = pl.col("filepath").map_elements(lambda s: os.path.basename(s)))
        .join(conversion, on="filename")
        .with_columns(new_filepath = pl.col("SampleID__").map_elements(lambda s: os.path.join(args.output_dir, s + ".fq.gz")))
        .select("filepath", "new_filepath")
    )
    logging.info(f"Found {samples.height} samples in {args.input_dir}")

    for row in samples.iter_rows():
        logging.info(f"Copying {row[0]} to {row[1]}")
        shutil.copyfile(row[0], row[1])

    logging.info("Done")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
