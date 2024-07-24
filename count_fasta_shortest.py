###############################
### count_fasta_shortest.py ###
###############################
# Author: Samuel Aroney
# Find the shortest contig
# python count_fasta_shortest.py [INPUT]

import sys
from Bio import SeqIO

input = sys.argv[1]

shortest = 10**5
with open(input) as input_file:
    for record in SeqIO.parse(input_file, 'fasta'):
        if len(record.seq) < shortest:
            shortest = len(record.seq)

print(input + "\t" + str(shortest))
