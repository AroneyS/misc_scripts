###########################
### find_exact_match.py ###
###########################
# Author: Samuel Aroney
# Find exact match sequence within fasta
# python find_exact_match.py [SEARCH SPACE FASTA] [SEED SEQUENCE FASTA] [OUTPUT]

import sys
from Bio import SeqIO

haystacks_path = sys.argv[1]
needles_path = sys.argv[2]
output_path = sys.argv[3]
delim = "\t"

haystacks = list(SeqIO.parse(haystacks_path, "fasta"))

matches = []
with open(needles_path) as needles_file:
    for needle in SeqIO.parse(needles_file, "fasta"):
        for haystack in haystacks:
            if needle.seq == haystack.seq:
                matches.append((needle.id, haystack.id))

if len(matches) > 0:
    with open(output_path, "w") as output_file:
        output_file.write("needle_id" + delim + "haystack_id" + "\n")
        for match in matches:
            output_file.write(match[0] + delim + match[1] + "\n")

print(f"Completed search for {needles_path} in {haystacks_path}")
