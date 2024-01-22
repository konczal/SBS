import sys
from Bio import SeqIO

fh = open(sys.argv[1], "r")

for seq_record in SeqIO.parse(fh, "fasta"):
     species_name = seq_record.id
     a = open(f"{species_name}.fa", "w")
     a.write(seq_record.format("fasta"))
     a.close()

fh.close()
