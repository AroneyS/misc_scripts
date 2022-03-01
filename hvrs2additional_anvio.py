################################
### hvrs2additional_anvio.py ###
################################
# Author: Samuel Aroney

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="find_HVRs output + splits_basic_info to additional_layers.")
parser.add_argument("--hvrs-output", type=str, metavar="<TSV FILE>", help="path to find_HVRs hvrs.tsv output file", required=True)
parser.add_argument("--splits-info", type=str, metavar="<TSV FILE>", help="path to Anvi'o splits_basic_info file", required=True)
parser.add_argument("--output", type=str, metavar="<TSV OUTPUT>", help="path to output file", required=True)

args = parser.parse_args()
hvrs_path = getattr(args, "hvrs_output")
splits_path = getattr(args, "splits_info")
output_path = getattr(args, "output")

hvrs_output = pd.read_csv(hvrs_path, sep = "\t")
hvrs_output["hvr_mid"] = hvrs_output.apply(lambda row: (row['hvr_start'] + row['hvr_end']) / 2, axis = 1)

splits = pd.read_csv(splits_path, sep = "\t")

merged = pd.merge(splits, hvrs_output, left_on = "parent", right_on = "contig")
output = merged[(merged["hvr_mid"] - merged["start"] >= 0) & (merged["end"] - merged["hvr_mid"] >= 0)]
output = output.groupby(["split", "sample"], as_index = False)["sample"].count().groupby(["split"], as_index = False)["sample"].count()
output.rename(columns = {"sample":"HVR_samples"}, inplace = True)
output.to_csv(output_path, sep="\t", index=False, quoting=0)
