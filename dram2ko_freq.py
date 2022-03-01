#######################
### dram2ko_freq.py ###
#######################
# Author: Samuel Aroney
# Convert DRAM annotations to KO frequency matrix

import argparse, os
import pandas as pd

class DramAnnotationConverter:
    ANNOTATIONS_FILENAME = "annotations.tsv"
    def __init__(self, annotations_dir, output_path):
        annotations_files = [os.path.join(dir, file) for dir,_,files in os.walk(annotations_dir, followlinks=True) for file in files if file == "annotations.tsv"]
        annotations = [self._load_annotation(file) for file in annotations_files]
        self.annotations = self._combine_and_clean(annotations)

        self.output_path = output_path
    
    def _load_annotation(self, filepath):
        annotation_df = pd.read_csv(filepath, sep = "\t", index_col = 0)
        return annotation_df
    
    def _combine_and_clean(self, annotations):
        df = pd.concat(annotations)
        # remove secondary ID annotations
        df["kegg_id"].replace(',.*','', regex=True, inplace=True)
        return df
    
    def convert(self):
        kegg_annotations = self.annotations.value_counts(subset=["fasta", "kegg_id"], dropna=True).to_frame().reset_index()
        kegg_frequency_matrix = kegg_annotations.pivot_table(index="kegg_id", columns="fasta", fill_value=0).reset_index()
        kegg_frequency_matrix.columns = kegg_frequency_matrix.columns.droplevel()
        self.kegg_frequency_matrix = kegg_frequency_matrix.rename(columns={"":"ID"})
    
    def output(self):
        self.kegg_frequency_matrix.to_csv(self.output_path, sep = "\t", index = False)

def main():
    parser = argparse.ArgumentParser(description="Convert DRAM annotations to KO frequency matrix.")
    parser.add_argument("--annotations", type=str, metavar="<DIR>", help="path to DRAM genome annotations directory", required=True)
    parser.add_argument("--output", type=str, metavar="<OUTPUT>", help="path to output file (tsv)", required=True)

    args = parser.parse_args()
    annotations_dir = getattr(args, "annotations")
    output_path = getattr(args, "output")

    converter = DramAnnotationConverter(annotations_dir, output_path)
    converter.convert()
    converter.output()

    
if __name__ == "__main__":
    main()
