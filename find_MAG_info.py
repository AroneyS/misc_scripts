##############################
### filter_fasta_length.py ###
##############################
# Author: Samuel Aroney
# Find MAG size, N-contigs and N-contigs >= 3kb
# python filter_fasta_length.py [INPUT]

import sys
import os
from Bio import SeqIO

input = sys.argv[1]

size = 0
n_contigs = 0
n_contigs_3kb = 0

with open(input) as input_file:
    for record in SeqIO.parse(input_file, 'fasta'):
        size += len(record.seq)
        n_contigs += 1
        if len(record.seq) >= 3000:
            n_contigs_3kb += 1

mag_name = os.path.basename(input).split(".")[0]
print(f"{mag_name}\t{size}\t{n_contigs}\t{n_contigs_3kb}")
