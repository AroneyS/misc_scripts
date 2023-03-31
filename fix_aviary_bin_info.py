#!/usr/bin/env python3

"""
Author: Samuel Aroney
Create bin_info file for Aviary recover without coverage
"""

import os
import sys
import argparse
import logging
import pandas as pd

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument("--checkm1", type=argparse.FileType("r"), help="CheckM1 results file")
    parser.add_argument("--checkm2", type=argparse.FileType("r"), help="CheckM2 results file", required=True)
    parser.add_argument("--bac-summary", type=argparse.FileType("r"), help="GTDBtk bacteria file", required=True)
    parser.add_argument("--ar-summary", type=argparse.FileType("r"), help="GTDBtk archaea file", required=True)
    parser.add_argument("--output", type=argparse.FileType("w"), help="Output bin_stats file", required=True)
    parser.add_argument("--output-minimal", type=argparse.FileType("w"), help="Output checkm_minimal file")

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

    bac_summary = pd.read_csv(args.bac_summary, sep="\t")
    ar_summary = pd.read_csv(args.ar_summary, sep="\t")
    taxa = pd.concat([bac_summary, ar_summary])
    taxa.rename({"user_genome" : "Bin Id"}, inplace=True, axis=1)

    checkm2 = pd.read_csv(args.checkm2, sep="\t")
    checkm2.rename({"Name" : "Bin Id", "Completeness" : "Completeness (CheckM2)", "Contamination" : "Contamination (CheckM2)"}, inplace=True, axis=1)

    if args.checkm1:
        checkm1 = pd.read_csv(args.checkm1, sep="\t", comment="[")
        checkm1.rename({"Completeness" : "Completeness (CheckM1)", "Contamination" : "Contamination (CheckM1)"}, inplace=True, axis=1)
        checkm_both = pd.merge(checkm1, checkm2, on=[checkm1.columns[0]])
    else:
        checkm_both = checkm2

    merged_out = pd.merge(checkm_both, taxa, on=[checkm_both.columns[0]])
    merged_out.to_csv(args.output, sep='\t', index=False)

    if args.output_minimal:
        checkm_minimal = checkm_both[["Bin Id",  "Marker lineage",  "# genomes", "# markers", "# marker sets",
                                        "0", "1", "2", "3", "4", "5+", "Completeness (CheckM1)", "Contamination (CheckM1)",
                                        "Completeness (CheckM2)", "Contamination (CheckM2)", "Strain heterogeneity"]]

        checkm_minimal.to_csv(args.output_minimal, sep="\t", index=False)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
