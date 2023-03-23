#!/usr/bin/env python
# Author: Samuel Aroney
# Remove duplicate reads from a FASTQ file
# python deduplicate_fastq.py [INPUT FASTQ] [OUTPUT FASTQ]

import sys
from Bio import SeqIO

input_path = sys.argv[1]
output_path = sys.argv[2]

detected = set()
unique = []
for rec in SeqIO.parse(open(input_path), "fastq"):
    if rec.description not in detected:
        detected.add(rec.description)
        unique.append(rec)

SeqIO.write(unique, open(output_path, "w"), "fastq")
