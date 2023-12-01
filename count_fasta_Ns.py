##############################
### count_fasta_Ns.py ###
##############################
# Author: Samuel Aroney
# Count Ns per contig
# python count_fasta_Ns.py [INPUT]

import sys
from Bio import SeqIO

input = sys.argv[1]

with open(input) as input_file:
    for record in SeqIO.parse(input_file, 'fasta'):
        print(record.id + "\t" + str(record.seq.count('N')))
