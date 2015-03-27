#!/usr/bin/python

import re
import sys
import os
import argparse
import sqlite3
from Bio import SeqIO
from Bio.Seq import Seq
import pysam 
from collections import Counter, OrderedDict


from utils.ucparse import uclustUserParser
from utils.cigar import cleanupCigar
from utils.samformater import makeSAMHdr, generateReadGroups
from utils.samToBam import samToBam
from utils.sqlitedb import buildsqlitedb
from utils.consensus import updateConsensus
from utils.utils import removeFiles, readfq, CONSENSUS_NAME
from utils.threadpool import ProducerConsumer

SAMHDR1 = """@HD\tVN:1.0"""
SAMHDR2 = """@SQ\tSN:%(rname)s\tLN:%(rlen)s"""
READGROUP = """@RG\tID:%(groupid)s\tSM:%(library)s\tLB:%(library)s\tPL:ILLUMINA"""

SAMLINE_S="""%(qname)s\t0\t%(rname)s\t%(pos)s\t0\t%(cigar)s\t*\t0\t0\t%(seq)s\t%(qual)s\tRG:Z:%(rgroup)s"""
SAMLINE_M="""%(qname)s\t0\t%(rname)s\t%(pos)s\t0\t%(cigar)s\t*\t0\t0\t%(seq)s\t%(qual)s"""

def makeSAMrecSingle(pos, seqinfo, cig, orient, refname, single = False, grpident = None):
    newcig, pos, rmfront, rmback = cleanupCigar(pos, cig, len(seqinfo[1]))
    
    if orient == '+':
        seq= seqinfo[1]
        qual = seqinfo[2]
    else:
        seq = str(Seq(seqinfo[1]).reverse_complement())          
        qual = seqinfo[2][::-1]

    if rmback:
        seq = seq[:-rmback]
        qual = qual[:-rmback]
    if rmfront:
        seq = seq[rmfront:]
        qual = qual[rmfront:]
    pos =  pos + 1 # need to shift from 0 base to 1 base number system
    if single:
        return SAMLINE_S%dict(qname = seqinfo[0], rname = CONSENSUS_NAME, pos = pos, cigar = newcig, seq = seq, qual = qual, rgroup = "GROUP-%s"%(grpident))
    else:
        return SAMLINE_M%dict(qname = seqinfo[0], rname = CONSENSUS_NAME, pos = pos, cigar = newcig, seq = seq, qual = qual)



def makeSAMrec(pos, seqinfo, cig, orient, refname, single = False, grpident = None):
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
        if grpident == None:
            tags.append( ("RG", "GROUP-%s"%(seqinfo.id.split("_")[0]),) )
        else:
            tags.append( ("RG", "GROUP-%s"%(grpident),) )
        a.tags = tags
    return a


def singleProducer(args):
    conseq = args[0]    
    clusterdat = args[1]
    idx = args[2]
    r = args[3]
    srec = []
    grps = dict()
    grpHdrs = dict()
    readgroups = []
    gidx = 0
    for d in clusterdat:
        grpident = d[0].split("_")[0]
        if grpident not in grps:
            grps[ grpident ] = gidx
            gid = "GROUP-%s"%(gidx)
            grpHdrs[grpident] = READGROUP%dict(groupid = gid, library = grpident )
            readgroups.append(gid)
            gidx += 1
        srec.append(makeSAMrecSingle(0, d[3], d[1], d[2], conseq[0], True, grps[ grpident ]))
    fname = "outFile_%s_%s"%(idx, len(conseq[1]) )
    tmp = str(conseq[1])   
    #rgrps, idtomap = generateReadGroups(grps.iterkeys())
    samdat ="%s\n%s\n%s\n%s\n"%(SAMHDR1, SAMHDR2%dict(rname = CONSENSUS_NAME, rlen = len(conseq[1]) ), "\n".join(grpHdrs.itervalues()), "\n".join(srec) )
    bamname, bamidxname = samToBam(samdat, fname, False)
    consensus = updateConsensus(bamname)
    bamdat = open(bamname).read() 
    bamidx = open(bamidxname).read() 
    removeFiles([bamname, bamidxname])    
    return (fname, readgroups, samdat, bamdat, bamidx, consensus,)

    
def singleConsumer(con, returndata):
    cur = con.cursor()
    for data in returndata:
        fname = data[0]
        rgrps = data[1]
        samdat = data[2]
        bamdat = data[3]
        bamidx = data[4]
        consensus = data[5]
        cur.execute("INSERT INTO files(name) values(?)""", (fname,) )
        cur.execute("""SELECT id FROM files WHERE name = ?;""", (fname,))
        row = cur.fetchone()
        fileID = row[0]
        cur.executemany("""INSERT INTO groups(fileID, groupid) VALUES(?,?);""", ( (fileID, g,) for g in  rgrps)  )      
        #cur.executemany("""INSERT INTO groups(fileID, groupid) VALUES(?,?);""", ( (fileID, g['ID'],) for g in  rgrps)  )      
        cur.execute("""INSERT INTO trimmed_inferSAM(fileID, sam, bam, bamidx) VALUES(?,?,?,?);""", (fileID, samdat, sqlite3.Binary(bamdat), sqlite3.Binary(bamidx),) )
        cur.execute("""INSERT INTO trimmed_consensus(fileID, sequence) VALUES(?, ?)""", (fileID, consensus,) )
    con.commit()
    


def processSingle(args):
    con = buildsqlitedb(args.database, False)
    cur = con.cursor()
    parser = uclustUserParser(args.ucinput)
    cleanids = set([ i.strip() for i in open(args.clean)])
    seqindex =  dict( ( (name, (name,seq,qual,), )    for name, seq, qual in readfq(open(args.sequences)) ) )
    conindex =  dict( ( (name, (name,seq, ), )    for name, seq, qual in readfq(open(args.consensus)) if (name in cleanids and args.minlen <= len(seq)) ) )
    cleanids = set(conindex.iterkeys())
    cutoff= int(args.min_clust)
    maxcutoff = int(args.max_clust)
    parser.parse(seqindex, cleanids, cutoff, maxcutoff)
    clusterid = parser.getClusters().keys()
    clusterid.sort()
    #for idx, r in enumerate(clusterid):
        #singleProducer((conindex[r], parser.getClusters(r), idx, r,))        
    worker = ProducerConsumer(args, args.threads, singleProducer, singleConsumer)
    worker.run( con, ( (conindex[r], parser.getClusters(r), idx, r,)  for idx, r in enumerate(clusterid) ) )
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
    removeFiles([samname])    


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
    parser_sub.add_argument('-t', '--threads', required = False, default = 1, type = int, help = "Processing threads(default: 1)")
    parser_sub.add_argument('-s', '--minlen', required = False, default = 0, type = int, help = "Minimum Sequence length(default: 0)")
    parser_sub.set_defaults(func = processSingle)
    args = parser.parse_args()
    
    args.func(args)
