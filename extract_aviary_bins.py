#!/usr/bin/env python3

"""
Author: Samuel Aroney
Extra quality bins from Aviary output, rename and collate metadata
"""

import sys
import os
import argparse
import logging
import extern
import polars as pl


def copy_file(source, dest):
    logging.info(f"Copying {source} to {dest}")

    os.makedirs(os.path.dirname(dest), exist_ok=True)
    cmd = f"cp {source} {dest}"
    extern.run(cmd)

    return True

def pipeline(aviary_dir, output_dir, prefix, min_completeness, max_contamination):
    bins = pl.read_csv(f"{aviary_dir}/bins/bin_info.tsv", separator = "\t"
        ).select(
            pl.col("Bin Id").alias("bin_id"),
            pl.col("Completeness (CheckM2)").alias("completeness"),
            pl.col("Contamination (CheckM2)").alias("contamination"),
            "Completeness_Model_Used", "classification", "classification_method",
            "translation_table", "red_value"
        ).filter(
            (pl.col("completeness") >= min_completeness) & (pl.col("contamination") <= max_contamination)
        ).with_row_count().with_columns(
            bin_name = pl.lit(prefix) + pl.lit("_") + pl.col("row_nr").cast(str),
        ).with_columns(
            bin_source = pl.col("bin_id").apply(lambda x: os.path.join(f"{aviary_dir}/bins/final_bins/", x + ".fna")),
            bin_dest = pl.col("bin_name").apply(lambda x: os.path.join(output_dir, x + ".fna")),
        ).with_columns(
            pl.struct(["bin_source", "bin_dest"]).apply(lambda x: copy_file(x["bin_source"], x["bin_dest"])).alias("copied")
        ).drop([
            "row_nr", "bin_source", "bin_dest", "copied"
        ])

    return bins

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument('--aviary-recover', help='Aviary recover directory', required=True)
    parser.add_argument("--output", help="Output folder", required=True)
    parser.add_argument("--prefix", help="Bin naming prefix", required=True)
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

    output = pipeline(args.aviary_recover, args.output, args.prefix, args.min_completeness, args.max_contamination)
    metadata_path = os.path.join(args.output, args.prefix + "_metadata.tsv")
    logging.info(f"Writing metadata to {metadata_path}")
    output.write_csv(metadata_path, separator = "\t")

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
