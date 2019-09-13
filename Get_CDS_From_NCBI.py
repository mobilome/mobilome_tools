import time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez

Entrez.email = "50883553@qq.com"

accession = []
with open("TIGD/TIGD5.fasta") as f:
    file = f.readlines()
for line in file:
    if line.startswith(">"):
        id = line.split("|")[3]
        accession.append(id)
for id in accession:
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
    records = SeqIO.read( handle, 'genbank')
    handle.close()
    features = records.features
    for feature in features:
        if feature.type == "CDS":
            rec_CDS =  SeqRecord(feature.location.extract(records).seq)
            rec_CDS.id = records.id
            rec_CDS.description = "CDS"
    with open("TIGD5_CDS.fasta", "a") as output_handle:
        SeqIO.write(rec_CDS, output_handle, "fasta")
    time.sleep(0.5)
