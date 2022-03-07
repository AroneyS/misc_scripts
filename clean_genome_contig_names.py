####################################
### clean_genome_contig_names.py ###
####################################
# Author: Samuel Aroney
# Remove underscores from genome and contig names
# Underscore separate genome and contig names
# python clean_genome_contig_names.py [GENOME_NAME] [INPUT] [OUTPUT]

import sys
from Bio import SeqIO
from BCBio import GFF

genome = sys.argv[1]
input = sys.argv[2]
output = sys.argv[3]

print(f"Cleaning genome {genome} contig names")

if input.endswith(".gff"):
    print(f"Cleaning GFF file")
    with open(output, "w") as output_file:
        with open(input) as input_file:
            for record in GFF.parse(input_file, target_lines=1):
                clean_record = genome + "_" + record.id.replace(genome + "_", "").replace("_", "-")
                print(f"Record {record.id} cleaned up to form {clean_record}")

                record.id = clean_record
                GFF.write([record], output_file)
                import pdb
                pdb.set_trace()

else:
    print(f"Cleaning fasta file")
    with open(output, "w") as output_file:
        with open(input) as input_file:
            for record in SeqIO.parse(input_file, 'fasta'):
                clean_record = genome + "_" + record.id.replace(genome + "_", "").replace("_", "-")
                print(f"Record {record.id} cleaned up to form {clean_record}")

                record.id = clean_record
                record.description = clean_record
                SeqIO.write(record, output_file, 'fasta')

