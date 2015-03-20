#!/usr/bin/python

import re
import sys
import os
import argparse
import sqlite3
from Bio import SeqIO
from Bio.Seq import Seq
import pysam 
from collections import Counter

from utils.ucparser import uclustUserParser
from utils.cigar import cleanupCigar
from utils.samformater import makeSAMHdr, generateReadGroups
from utils.sqlitedb import buildsqlitedb
from utils.consensus import updateConsensus
from utils.utils import removeFiles


def makeSAMrec(pos, seqinfo, cig, orient, refname, single = False):
    newcig, pos, rmfront, rmback = cleanupCigar(pos, cig, len(seqinfo.seq))
    a = pysam.AlignedRead()
    a.tid = 0
    a.rname = 0
    a.qname = seqinfo.id   
    a.flag = 0x00 
    if orient == '-':
        a.seq = str(seqinfo.seq.reverse_complement())
    else:
        a.seq= str(seqinfo.seq)
    if rmback:
        a.seq = a.seq[:-rmback]
    if rmfront:
        a.seq = a.seq[rmfront:]

    if seqinfo.letter_annotations:
        tmpq = seqinfo.format("fastq").strip().split("\n")[-1]
        if orient == '-':
            tmpq = tmpq[::-1]
        if rmback:
            tmpq = tmpq[:-rmback]
        if rmfront:
            tmpq = tmpq[rmfront:]
        a.qual = tmpq

    a.pos =  pos
    a.cigarstring = newcig   
    a.rnext = -1
    a.pnext= -1
    a.tlen = 0
    if single:
        tags=  []
        tags.append( ("RG", "GROUP-%s"%(seqinfo.id.split("_")[0]),) )
        a.tags = tags
    return a


def processSingle(args):
    con = buildsqlitedb(args, args.database, False)
    cur = con.cursor()
    parser = uclustUserParser(args.ucinput)
    seqindex = SeqIO.index(args.sequences, "fastq")
    conindex = SeqIO.index(args.consensus, "fasta")
    cutoff= int(args.min_clust)
    maxcutoff = int(args.max_clust)
    parser.parse(seqindex, set([ i.strip() for i in open(args.clean)]), cutoff, maxcutoff)
    srec = []
    reflen = 0
    start = 0
    rseq = ""
    clusterid = parser.getClusters().keys()
    clusterid.sort()
    for idx, r in enumerate(clusterid):
        conseq = conindex[r]
        
        grps = set()
        for d in parser.getClusters()[r]:
            srec.append(makeSAMrec(start, d[3], d[1], d[2], conseq.id, True))
            grps.add(d[0].split("_")[0])
            if d[-1]:
                tmp = str(d[3].seq)
                o.write(tmp)
                print >> oo, ">%s\n%s"%(d[3].id, tmp)

        fname = "outFile_%s_%s"%(idx, len(conseq.seq) )
        cur.execute("INSERT INTO files(name) values(?)""", (fname,) )
        cur.execute("""SELECT id FROM files WHERE name = ?;""", (fname,))
        row = cur.fetchone()
        fileID = row[0]
        tmp = str(conseq.seq)
        
        rgrps, idtomap = generateReadGroups(grps)
        cur.executemany("""INSERT INTO groups(fileID, groupid) VALUES(?,?);""", ( (fileID, g['ID'],) for g in  rgrps)  )      
        samout = pysam.Samfile("%s.sam"%(fname ) , "wh", header= makeSAMHdr("Consensus", len(conseq.seq), rgrps))
        for f in srec:
            samout.write(f)
        samout.close()
        samname = "%s.sam"%(fname)
       
        samdat = open(samname).read()
        bamname, bamidxname = samToBam(samdat, prefix, False)
        tmp = updateConsensus(bamname)
        bamdat = open(bamname).read() 
        bamidx = open(bamidxname).read() 
        removeFiles([samname, bamname, bamidxname])

        cur.execute("""INSERT INTO trimmed_inferSAM(fileID, sam, bam, bamidx) VALUES(?,?,?,?);""", (fileID, samdat, sqlite3.Binary(bamdat), sqlite3.Binary(bamidx),) )
        cur.execute("""INSERT INTO trimmed_consensus(fileID, sequence) VALUES(?, ?)""", (fileID, tmp,) )

        srec = []
    con.commit()
    con.close()


def processMultiple(args):
    parser = uclustUserParser(args.ucinput)
    seqindex = SeqIO.index(args.sequences, "fastq")
    conindex = SeqIO.index(args.consensus, "fasta")
    cutoff = int(args.min_clust)
    maxcutoff = int(args.max_clust)
    prefix = args.outprefix
    refname = prefix
    parser.parse(seqindex, set([ i.strip() for i in open(args.clean)]), cutoff, maxcutoff)
    srec = []
    reflen = 0
    start = 0
    rseq = ""
    clusterid = parser.getClusters().keys()
    clusterid.sort()

    oo = open("%s_pseudo_ref_parts.fasta"%(prefix), "w")
    with open("%s_pseudo_ref.fasta"%(prefix), "w") as o:
        print >> o, ">%s"%(refname)
        for r in clusterid:
            start = reflen
            for d in parser.getClusters()[r]:
                srec.append(makeSAMrec(start, d[3], d[1], d[2], refname))
                if d[-1]:
                    tmp = str(d[3].seq)
                    o.write(tmp)
                    print >> oo, ">%s\n%s"%(d[3].id, tmp)
                    reflen += len(tmp) 
            conseq = conindex[r]
            tmp = str(conseq.seq)
            o.write(tmp)
            print >> oo, ">%s\n%s"%(conseq.id, tmp)
            reflen += len(tmp) 
    oo.close()                
    samout = pysam.Samfile("%s.sam"%(prefix) , "wh", header= makeSAMHdr(refname, reflen))
    for f in srec:
        samout.write(f)
    samout.close()
    samname = "%s.sam"%(prefix)
  
    bamname, bamidxname = samToBam(open(samname).read(), prefix, False)

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-u', '--ucinput', required = True, help = "userout generated by vsearch/usearch" )
    parser.add_argument('-q', '--sequences',  required = True, help = "fastq file of input sequences")
    parser.add_argument('-c', '--clean', required = True, help = "Clean ID list file")
    parser.add_argument('-f', '--consensus',  required = True, help = "Fasta format of the consensus sequences")
    parser.add_argument('-l', '--min_clust', required = False, default = 1, type = int, help = "Minimum cluster size (default: 1)")
    parser.add_argument('-m', '--max_clust', required = False, default = 0, type = int,  help = "Maximum cluster size (default: infinite)")

    subparsers = parser.add_subparsers(dest='action', help='Available commands')

    parser_sub = subparsers.add_parser('multiple')
    parser_sub.add_argument('-o', '--outprefix', required = True, help = "Prefix for all output files")
    parser_sub.set_defaults(func = processMultiple)

    parser_sub = subparsers.add_parser('single')
    parser_sub.add_argument('-d', '--database', required = True, help = "Database name")
    parser_sub.set_defaults(func = processSingle)
    args = parser.parse_args()
    
    args.func(args)
