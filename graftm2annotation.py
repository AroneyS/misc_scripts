############################
### graftm2annotation.py ###
############################
# Author: Samuel Aroney
# GraftM *_read_tax.tsv to iTOL annotation file

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="GraftM *_read_tax.tsv to iTOL annotation file.")
parser.add_argument("--graftm-output", type=str, metavar="<TSV FILE>", help="path to graftm *_read_tax.tsv output file", required=True)
parser.add_argument("--grouping", type=str, metavar="<TSV FILE>", help="path to grouping and colour file", required=True)
parser.add_argument("--output", type=str, metavar="<TSV OUTPUT>", help="path to output file", required=True)

args = parser.parse_args()
graftm_output_path = getattr(args, "graftm_output")
grouping_path = getattr(args, "grouping")
output_path = getattr(args, "output")

graftm_output = pd.read_csv(graftm_output_path, sep="\t", header=None, names = ["id", "type"])
grouping = pd.read_csv(grouping_path, sep="\t")

output = pd.merge(graftm_output, grouping, on = "type", how = "left")
output["id"] =  output["id"].apply(lambda x: "'" + str(x) + "'")
output.loc[:, ["id", "colour", "group"]].to_csv(output_path, sep="\t", index=False, header=False, quoting=0)
