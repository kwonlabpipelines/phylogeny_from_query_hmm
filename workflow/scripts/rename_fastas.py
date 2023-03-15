import os
from Bio import SeqIO

new_recs = []
with open(snakemake.input[0], "r") as f:
    species = snakemake.wildcards.genome
    recs = list(SeqIO.parse(f, "fasta"))
    for i,rec in enumerate(recs):
        rec.id = "{}_{:03d}".format(species, i)
        rec.description = ""
        new_recs.append(rec)

with open(snakemake.output[0], "w") as f:
    SeqIO.write(new_recs, f, "fasta")