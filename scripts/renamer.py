#!/usr/bin/python

from Bio import SeqIO
import sys

with open(sys.argv[3], "w") as o:
    for s in SeqIO.parse(sys.argv[1], "fastq"):
        s.id = "%s_%s"%(sys.argv[2], s.id)
        SeqIO.write([s], o, "fastq")
