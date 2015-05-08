#!/usr/bin/env python2.7


#arr=("l24" "plob1" "plob2" "plob3" "pcom1" "pcom2" "pcom3");for i in ${arr[@]}; do echo ${i}; findOrder.sh ${i}_pseudo_ref_parts.fasta ${i} 2 21; done
#findBestOrderv2.py . _kmer.filter order.lst 


from Bio import SeqIO
import sys
from collections import defaultdict

if len(sys.argv) != 4:
    print "USAGE: %s <kmer counts> <output prefix> <min identity>"%(sys.argv[0])
    sys.exit(1)

kcounts = defaultdict(int)
minident = int(sys.argv[3])

total_kmers = 0
output = "%s_kmer.filter"%(sys.argv[2])
with open("%s.tmp"%(output), "w") as o:
    infile = open(sys.argv[1])
    for s in infile:
        ident = int(s[1:])
        if ident < minident:
            infile.next()
            continue
        total_kmers += ident
        kcounts[ident] += 1
        seq = infile.next().strip()
        print >> o,  "%s\t%s"%(seq, ident)

total_kmers = float(total_kmers)
discards = set(xrange(minident))

with open(output, "w") as o:
    for s in open("%s.tmp"%(output)):
        s = s.strip().split("\t")
        ident = int(s[1])
        if ident in discards:
            continue
        print >> o,  "%s\t%s"%(s[0], (float(ident) / total_kmers))

try:
    os.remove("%s.tmp"%(output))
except:
    pass

