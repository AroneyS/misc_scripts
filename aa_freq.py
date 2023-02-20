##################
### aa_freq.py ###
##################
# Author: Samuel Aroney
# Count amino acids across fasta file
# python aa_freq.py [FASTA]

import sys
import os
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

input = sys.argv[1]
aminoacids = defaultdict(int)

with open(input) as f:
    for record in SeqIO.parse(f, 'fasta'):
        for aa in record.seq:
            aminoacids[aa] +=1

aminoacids_df = pd.DataFrame([(a, aminoacids[a]) for a in aminoacids], columns = ["aminoacid", "count"])
aminoacids_df["freq"] = aminoacids_df["count"] / sum(aminoacids_df["count"])
aminoacids_df["filename"] = os.path.basename(input)

aminoacids_df.to_csv(sys.stdout, columns=["filename", "aminoacid", "freq"], header=False, index=False, sep="\t")
