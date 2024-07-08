#!/usr/bin/env python3

"""
Author: Samuel Aroney
Extract quality bins from Bin Chicken output, rename and collate metadata
"""

import sys
import os
import argparse
import logging
import extern
import polars as pl


def copy_file(source, dest):
    logging.debug(f"Copying {source} to {dest}")

    os.makedirs(os.path.dirname(dest), exist_ok=True)
    cmd = f"cp {source} {dest}"
    extern.run(cmd)

    return True

def pipeline(binchicken_dir, output_dir, single_sample=False):
    clusters = (
        pl.read_csv(f"{binchicken_dir}/coassemble/target/elusive_clusters.tsv", separator="\t")
        .select("coassembly", "samples")
    )

    if single_sample:
        clusters = (
            clusters
            .with_columns(
                prefix = pl.col("samples")
                    .str.replace("MainAutochamber", "MA")
                    .str.replace("PalsaHole", "PH")
                    .str.replace("Hodgkins", "H")
                    .str.replace("TempControlInc", "TCI")
                    .str.replace(r"^", "cmr1.")
                )
        )
    else:
        logging.error("Coassembly mode not implemented yet")
        raise NotImplementedError

    bins = (
        clusters
        .with_columns(
            bin_info = pl.col("coassembly").map_elements(
                lambda x: os.path.join(binchicken_dir, "coassemble", "coassemble", x, "recover", "bins", "bin_info.tsv")
                )
            )
        .with_columns(
            pl.col("bin_info").map_elements(
                lambda x: pl.read_csv(x, separator="\t")
                    .select(
                        pl.col("Bin Id").alias("bin_id"),
                        pl.col("Completeness (CheckM2)").alias("completeness"),
                        pl.col("Contamination (CheckM2)").alias("contamination"),
                        "Completeness_Model_Used",
                        )
                    .to_struct("bin_info")
                )
            )
        .explode("bin_info")
        .unnest("bin_info")
        .with_row_count()
        .with_columns((pl.col("row_nr") - pl.min("row_nr")).over("coassembly"))
        .with_columns(bin_name = pl.concat_str(pl.col("prefix") + pl.lit(".") + pl.col("row_nr").cast(str)))
        .with_columns(
            bin_source = pl.struct(["coassembly", "bin_id"]).map_elements(
                lambda x: os.path.join(binchicken_dir, "coassemble", "coassemble", x["coassembly"], "recover", "bins", "final_bins", x["bin_id"] + ".fna")
                ),
            bin_dest = pl.col("bin_name").map_elements(lambda x: os.path.join(output_dir, "bins", x + ".fna")),
            )
        .with_columns(
            copied = pl.struct(["bin_source", "bin_dest"]).map_elements(lambda x: copy_file(x["bin_source"], x["bin_dest"]))
            )
        .drop(["row_nr", "bin_source", "bin_dest", "copied"])
    )

    return bins

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--debug", help="output debug information", action="store_true")
    parser.add_argument("--quiet", help="only output errors", action="store_true")

    parser.add_argument("--binchicken", help="Bin Chicken output directory", required=True)
    parser.add_argument("--output", help="Output folder", required=True)
    parser.add_argument("--single-sample", help="Single sample mode", action="store_true")

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format="%(asctime)s %(levelname)s: %(message)s", datefmt="%Y/%m/%d %I:%M:%S %p")

    output = pipeline(args.binchicken, args.output, args.single_sample)

    metadata_path = os.path.join(args.output, "metadata.tsv")
    logging.info(f"Writing metadata to {metadata_path}")
    output.write_csv(metadata_path, separator = "\t")

    logging.info("Done")

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
