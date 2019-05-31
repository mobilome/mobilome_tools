import sys
from Bio import SeqIO
input_name = "X:/seagate/pig/pig_genome/GC"
output_name = input_name + ".info"

f = open(input_name,'r')
out = open(output_name,'w')
for rec in SeqIO.parse(f, 'fasta'):
    out.write("%s\t%s\n" % (rec.id,len(rec.seq)))
f.close()
out.close()
