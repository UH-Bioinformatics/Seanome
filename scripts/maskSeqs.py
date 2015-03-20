#!/usr/bin/python                                                                                                                        
import argparse
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


def getLines(infile, qcov):
    lines = defaultdict(list)
    with open(infile) as f:
        for l in f:
            l = l.strip().split()
            if int(l[3]) < qcov:
                continue
            lines[l[1]].append(l)
    for k, v in lines.iteritems():
        v.sort(key = lambda x : int(x[4]) , reverse = True )
    return lines

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--id',  required = True, help = "bad ids output")
    parser.add_argument('-f', '--fasta', required = True, help = "masked contigs that make up the targets in the badids file")
    parser.add_argument('-c', '--coverage', type = float, default = 75, help = "Target Coverage threshold between drop and N mask (default 75) [0-100]" )
    parser.add_argument('-q', '--qcoverage', type = float, default = 75, help = "Query Coverage threshold.  Ignore hits that are below this threshold. (default 75) [0-100]" )
    parser.add_argument('-o', '--output', required = True, help = "output fasta file")
    args = parser.parse_args()
    
    lines = getLines(args.id, args.qcoverage)
    
    with open(args.output, "w") as o:
        for entry in SeqIO.parse(args.fasta, "fasta"):
            if entry.id not in lines:
                if (100.0 * (float(str(entry.seq).count("N")) / float(len(str(entry.seq))) )) >= args.coverage:
                    continue
                SeqIO.write([entry], o, "fasta")
            else:
                dat = lines[entry.id]
                #print dat
                if int(dat[0][4]) >= args.coverage:               
                    #print "drop"
                    continue
                else:
                    seq = str(entry.seq)
                    for r in dat:
                        start  = int(r[7]) - 1
                        end = int(r[8]) - 1
                        seq = "".join([ s if (i < start or i > end) else "N" for i, s in enumerate(seq)] )
                    entry.seq = Seq(seq)
                    if (100.0 * (float(seq.count("N")) / float(len(seq)))) >= args.coverage:
                        #print "More"
                        continue
                    #print seq
                    SeqIO.write([entry], o, "fasta")
