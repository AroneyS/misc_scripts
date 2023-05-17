#!/usr/bin/env python3

"""
Author: Samuel Aroney
Sum total cpu/wall time and max memory from provided PBS output files
"""

import sys
import argparse
import logging
import polars as pl
import re

def get_usage(pbs_out_path):
    with open(pbs_out_path, "r") as f:
        pbs_out = f.read()

    cpu_time = re.search("CPU time  : (.+)\n", pbs_out).group(1)
    wall_time = re.search("Wall time : (.+)\n", pbs_out).group(1)
    memory = re.search("Mem usage : (.+)\n", pbs_out).group(1)

    return cpu_time, wall_time, memory

def get_totals(pbs_out_paths):
    DURATION_REGEX=r"([+-]?\d+)"
    MEMORY_REGEX=r"(\d+)([kKmMgG]?[bB])"

    units_to_bytes = {
        "b": 1,
        "kb": 1024,
        "mb": 1024 ** 2,
        "gb": 1024 ** 3
    }

    df = (
        pl.DataFrame({"path": pbs_out_paths})
        .with_columns(
            data = pl.col("path").apply(get_usage)
            )
        .with_columns(
            pl.col("data").arr.to_struct(
                n_field_strategy="max_width",
                fields=["cpu_time", "wall_time", "mem_usage"]
                )
            )
        .unnest("data")
        .select(
            pl.col("cpu_time").str.extract_all(DURATION_REGEX),
            pl.col("wall_time").str.extract_all(DURATION_REGEX),
            pl.col("mem_usage").str.extract(MEMORY_REGEX).cast(float),
            pl.col("mem_usage").str.extract(MEMORY_REGEX, group_index=2).str.to_lowercase().alias("mem_unit")
            )
        .with_columns(
            pl.duration(
                hours=pl.col("cpu_time").arr.get(0),
                minutes=pl.col("cpu_time").arr.get(1),
                seconds=pl.col("cpu_time").arr.get(2)
                ).alias("cpu_time"),
            pl.duration(
                hours=pl.col("wall_time").arr.get(0),
                minutes=pl.col("wall_time").arr.get(1),
                seconds=pl.col("wall_time").arr.get(2)
                ).alias("wall_time"),
            (pl.col("mem_usage").cast(pl.Float64) * pl.col("mem_unit").apply(lambda unit: units_to_bytes[unit])).alias("mem_usage")
            )
        .select(
            pl.sum("cpu_time"),
            pl.sum("wall_time"),
            pl.max("mem_usage")
        )
    )

    return df

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument('--pbs-out', nargs='+', help='PBS output files', required=True)

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    output = get_totals(args.pbs_out)
    print(output)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
