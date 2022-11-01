#!/usr/bin/env python3

"""
Author: Samuel Aroney
Find info from Aviary jobs
"""

import sys
import argparse
import logging
from ruamel.yaml import YAML
import extern
import pandas as pd
import re


def load_config(config_file):
    yaml = YAML()
    yaml.version = (1, 1)
    yaml.default_flow_style = False

    return yaml.load(config_file)

def get_item(config, key):
    if key in config:
        return config[key]
    else:
        return None

def get_fq_size(path, cat = "zcat"):
    cmd = f"{cat} '{path}' | sed -n 2~4p | tr -d '\n' | wc -m"
    output = extern.run(cmd)

    return int(output.strip())

def get_usage(pbs_out_paths, output_dir, commands_file):
    """Load pbs result and match against commands list using output_dir
    """
    commands = pd.read_csv(commands_file, names = ["command"])
    matching_command = commands[commands["command"].str.contains(output_dir.strip("/"))]
    if len(matching_command) == 0:
        logging.error(f"Could not find matching command for output dir {output_dir}")
        sys.exit(1)
    elif len(matching_command) > 1:
        logging.error(f"Found multiple matching commands {str(matching_command.index.to_list())} for output dir {output_dir}")
        sys.exit(1)

    matching_pbs_out_path = sorted(pbs_out_paths)[matching_command.index[0]]
    with open(matching_pbs_out_path, "r") as f:
        pbs_out = f.read()

    cpu_time = re.search("CPU time  : (.+)\n", pbs_out).group(1)
    wall_time = re.search("Wall time : (.+)\n", pbs_out).group(1)
    memory = re.search("Mem usage : (.+)\n", pbs_out).group(1)

    return cpu_time, wall_time, memory

def get_info(output_dir, config_file, pbs_out_paths, commands_file):
    # Load config
    config = load_config(config_file)
    short_reads_1 = get_item(config, "short_reads_1")
    short_reads_2 = get_item(config, "short_reads_2")
    short_reads_size = 0
    if short_reads_1:
        for read in short_reads_1:
            short_reads_size += get_fq_size(read)
    if short_reads_2:
        for read in short_reads_2:
            short_reads_size += get_fq_size(read)

    long_reads = get_item(config, "long_reads")
    long_reads_size = 0
    if long_reads:
        for read in long_reads:
            long_reads_size += get_fq_size(read)

    assembly = get_item(config, "fasta")
    assembly_size = 0
    if assembly:
        assembly_size += get_fq_size(assembly[0], cat = "cat")

    # Load pbs result and match against commands list
    cpu_time, wall_time, memory = get_usage(pbs_out_paths, output_dir, commands_file)

    info = pd.DataFrame({
        "output_dir": [output_dir],
        "short_reads_size": [short_reads_size],
        "long_reads_size": [long_reads_size],
        "long_read_type": get_item(config, "long_read_type"),
        "assembly_size": [assembly_size] if assembly_size != 0 else [None],
        "cpu_time": [cpu_time],
        "wall_time": [wall_time],
        "memory": [memory]
    })

    return info

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument('--aviary-output', help='Aviary output directory', required=True)
    parser.add_argument('--config-file', type=argparse.FileType("r"), help='Aviary config file', required=True)
    parser.add_argument('--pbs-out', nargs='+', help='Aviary PBS output files', required=True)
    parser.add_argument('--commands', type=argparse.FileType("r"), help='Aviary commands file', required=True)

    parser.add_argument("-o", "--outfile", type=argparse.FileType("w"), default=sys.stdout, help="Output file (tsv)")
    parser.add_argument("--no-header", action="store_true", help="Output without header")

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    output = get_info(args.aviary_output, args.config_file, args.pbs_out, args.commands)

    output.to_csv(args.outfile, header = not args.no_header, index = False, sep = "\t")

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
