####################################
### clean_genome_contig_names.py ###
####################################
# Author: Samuel Aroney
# Remove underscores from genome and contig names
# Underscore separate genome and contig names
# python clean_genome_contig_names.py [GENOME_NAME] [INPUT] [OUTPUT]

import sys
from Bio import SeqIO
import gffutils

genome = sys.argv[1]
input = sys.argv[2]
output = sys.argv[3]

print(f"Cleaning genome {genome} contig names")

if input.endswith(".gff"):
    print(f"Cleaning GFF file")
    with open(output, "w") as output_file:
        db = gffutils.create_db(input, dbfn=output+".db", force=True, keep_order=True)
        output_file.write(f"## {db.directives[0]}\n")
        for record in db.all_features():
            clean_record = genome + "_" + record.id.replace(genome + "_", "").replace("_", "-")
            clean_contig = genome + "_" + record.seqid.replace(genome + "_", "").replace("_", "-")
            print(f"Record {record.id} from contig {record.seqid}, cleaned up to form {clean_record} from {clean_contig}")

            record.id = clean_record
            record.attributes["ID"] = clean_record
            record.attributes["locus_tag"] = clean_record
            record.seqid = clean_contig
            output_file.write(str(record) + "\n")

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

