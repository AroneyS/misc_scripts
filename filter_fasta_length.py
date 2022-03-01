##############################
### filter_fasta_length.py ###
##############################
# Author: Samuel Aroney
# Filter fasta file by length
# python filter_fasta_length.py [MIN LEN] [INPUT] [OUTPUT]

import sys
from Bio import SeqIO

n = int(float(sys.argv[1]))
input = sys.argv[2]
output = sys.argv[3]

with open(output, "w") as output_file:
    with open(input) as input_file:
        for record in SeqIO.parse(input_file, 'fasta'):
            if len(record.seq) >= n:
                SeqIO.write(record, output_file, 'fasta')

