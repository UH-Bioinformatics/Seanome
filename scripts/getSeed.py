#!/usr/bin/python

import sys
from Bio import SeqIO


dataQ = {}
dataT = {}
pos = 0
for seq in SeqIO.parse(sys.argv[2], "fasta"):
    end = pos + len(seq.seq)
    dataQ[seq.id] = pos
    pos = end

pos = 0
for seq in SeqIO.parse(sys.argv[3], "fasta"):
    end = pos + len(seq.seq)
    dataT[seq.id] =pos
    pos = end


seen = set()

with open(sys.argv[1]) as f:
    for l in f:
        l = l.strip().split()
        if l[0] in seen:
            continue
        seen.add(l[0])
        #query+target+id+alnlen+bits+qstrand+tstrand+qlo+qhi+tlo+thi+qrow+trow
        l[7]  = str(int(l[7]) + dataQ[l[0]])
        l[8]  = str(int(l[8]) + dataQ[l[0]])
        l[9]  = str(int(l[9]) + dataT[l[1]])
        l[10] = str(int(l[10]) + dataT[l[1]])
        print "\t".join(l)

