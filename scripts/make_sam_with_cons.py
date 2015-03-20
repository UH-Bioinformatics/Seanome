#!/usr/bin/python

import re
import sys
import os
import argparse
import sqlite3
from Bio import SeqIO
from Bio.Seq import Seq
import pysam 
from collections import defaultdict, Counter
import subprocess

class parser(object):
    def __init__(self, infile):
        self.infile = infile
        self.clusters = defaultdict(list)
        self.reps ={}


    def getClusters(self):
        return self.clusters


    def getReps(self):
        return self.reps


    def parse(self):
        pass


class uclustUserParser(parser):
    def __init__(self, infile):
        super(self.__class__, self).__init__(infile)


    def parse(self, index, use, cutoff = 0, maxcutoff = 0):
        if(not os.path.exists(self.infile)):
            print >> sys.stderr, "usearch did not generate an output.  Make sure usearch executed properly."
            return False
        with open(self.infile) as cfile:
            # sequence centroid cigar sequence_strand centroid_strand
            for l in cfile:
                l = l.strip().split()
                if l[1] not in use:
                    continue
                # a dedup sequence, just pull who it says it is a dup of..
                if l[5] == '*': 
                    self.clusters[l[1]].append( (l[0], l[2], l[3], index[l[0]], False,) )
                else:
                    self.clusters[l[1]].append( (l[0], l[2], l[3], index[l[5]], False,) )
        # strip out all clusters that do not meet the required cutoff
        for k in self.clusters.keys():
            if cutoff >= 0 and len(self.clusters[k]) < cutoff:
                self.clusters.pop(k, None)
            elif maxcutoff > 0 and len(self.clusters[k]) > maxcutoff:
                self.clusters.pop(k, None)
        return True


def interpretAln(parts):
    index = 0
    rm = 0
    pos = 0
    for idx, f in enumerate(parts):
        index = idx
        if f[1] == 'I':
            pos += int(f[0])
        elif f[1] == 'D':
            rm += int(f[0])
        else:
            break
    return index, rm, pos


FINDALL_CMP=re.compile(r'([0-9]*)([MID=])')
# I='D' or I='N'
CIGMAP = dict(M='M', D='I', I='D')

def cleanupCigar(pos, cig, seqlen):
    parts = FINDALL_CMP.findall(cig)
    # in all its wisdom.. usearch can return = for 
    # the alignment, if we have a perfect match.
    if len(parts) == 1:
        if parts[0][1] == '=':
            return "%sM"%(seqlen)
        else:
            return cig, pos, 0,0

    # fix the RLE so that everything has both elements
    parts = [ (p[0] if p[0] else '1', p[1],) for p in parts]

    f_idx, rmfront, shiftright = interpretAln(parts)
    # trim and reverse the parts so that we pull off the tail end for processing and dont overlap the forward
    r_idx, rmback, noshift = interpretAln(parts[f_idx:][::-1])

    # r_idx == 0 would make the final result empty
    if r_idx:
        parts = parts[f_idx: -r_idx]
    else:
        parts = parts[f_idx:]

    newcig = ''.join( ( "".join([p0, CIGMAP[p1]]) for p0,p1 in parts) )

    return newcig, pos + shiftright, rmfront, rmback 


def makeSAMrec(pos, seqinfo, cig, orient, refname, single = False):
    newcig, pos, rmfront, rmback = cleanupCigar(pos, cig, len(seqinfo.seq))
    a = pysam.AlignedRead()
    a.tid = 0
    a.rname = 0
    a.qname = seqinfo.id   
    #a.flag = 0x02 
    a.flag = 0x00 

    if orient == '-':
        #a.flag |= 0x10
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


def makeSAMHdr(ref, seqlen, rgroups = None):
    if rgroups:
        return dict(HD = dict(VN = '1.0'), SQ = [ {'SN': ref, 'LN': seqlen }], RG = rgroups )
    else:
        return dict(HD = dict(VN = '1.0'), SQ = [ {'SN': ref, 'LN': seqlen }])


def generateReadGroups(groups):
    RG_template = { 'ID': '',
                    'LB': '',
                    'SM': '',
                    'PL': 'ILLUMINA'} # Platform/technology used to produce the reads.  Spoof it                        
    readgroups = []
    idtomap = {}
    for idx, g in enumerate(groups):
        rgrp = RG_template.copy()
        rgrp['ID'] = "GROUP-%s"%(g) 
        rgrp['LB'] = g
        rgrp['SM'] = g
        readgroups.append(rgrp)
    return readgroups



def buildsqlitedb(args, dbname):
    con = sqlite3.connect(dbname)
    con.execute("""PRAGMA foreign_keys = ON;""")

    con.execute("""CREATE TABLE IF NOT EXISTS files(id INTEGER PRIMARY KEY, name TEXT);""")
    con.execute("""CREATE INDEX IF NOT EXISTS file_name_idx ON files(name ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS groups(id INTEGER PRIMARY KEY, fileID INTEGER, groupid TEXT, FOREIGN KEY(fileID) REFERENCES files(id) );""")
    con.execute("""CREATE INDEX IF NOT EXISTS groups_fileid_idx ON groups(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_consensus(id INTEGER PRIMARY KEY, fileID INTEGER, sequence TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_con_fileid_idx ON trimmed_consensus(fileID ASC);""")

    #con.execute("""CREATE TABLE IF NOT EXISTS inferSAM(id INTEGER PRIMARY KEY, fileID INTEGER, sam TEXT, bam BLOB, bamidx BLOB, FOREIGN KEY(fileID) REFERENCES files(id));""")
    #con.execute("""CREATE INDEX IF NOT EXISTS infersam_fileid_idx ON inferSAM(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_inferSAM(id INTEGER PRIMARY KEY, fileID INTEGER, sam TEXT, bam BLOB, bamidx BLOB, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_infersam_fileid_idx ON trimmed_inferSAM(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_vcf(id INTEGER PRIMARY KEY, fileID INTEGER, vcf TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_vcf_fileid_idx ON trimmed_vcf(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_modvcf(id INTEGER PRIMARY KEY, fileID INTEGER, vcf TEXT, json TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_modvcf_fileid_idx ON trimmed_modvcf(fileID ASC);""")

    con.commit()
    return con



def processSingle(args):
    con = buildsqlitedb(args, args.database)
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
        
        rgrps = generateReadGroups(grps)
        samout = pysam.Samfile("%s.sam"%(fname ) , "wh", header= makeSAMHdr("Consensus", len(conseq.seq), rgrps))
        for f in srec:
            samout.write(f)
        samout.close()
        samname = "%s.sam"%(fname)
        samdat, bamdat, bamidx, bamname, bamidxname = buildBAM(fname, samname)        
        cur.execute("""INSERT INTO trimmed_inferSAM(fileID, sam, bam, bamidx) VALUES(?,?,?,?);""", (fileID, samdat, sqlite3.Binary(bamdat), sqlite3.Binary(bamidx),) )
        cur.executemany("""INSERT INTO groups(fileID, groupid) VALUES(?,?);""", ( (fileID, g['ID'],) for g in  rgrps)  )
       
        tmp = updateConsensus(bamname)
        print tmp
        cur.execute("""INSERT INTO trimmed_consensus(fileID, sequence) VALUES(?, ?)""", (fileID, tmp,) )


        try:
            os.remove(bamname)
        except:
            pass
        try:
            os.remove(bamidxname)
        except:
            pass
        srec = []
    con.commit()
    con.close()



def updateConsensus(bamfile):
    samfile = pysam.Samfile(bamfile)
    colBases = []
    for pileupcolumn in samfile.pileup():
        bases =[]
        for pup in pileupcolumn.pileups:
            bases.append(pup.alignment.seq[pup.qpos])
        colBases.append(bases)
    newSeq=""
    for pos in colBases:
        values = sorted(Counter(pos).items(), key=lambda x: x[1], reverse=True)
        added = False
        for val in values:
            if val[0] == "N" or val[0] == "-":
                continue
            else:
                newSeq += val[0]
                added = True
                break
        if not added:
            newSeq += "N"
    return newSeq



def buildBAM(name, samout, store = True):
    samdat = open(samout, "rb").read()
    try:
        os.remove(samout)
    except:
        pass

    bamout = "%s.bam"%(name)
    bamidx = "%s.bai"%(bamout)
    cline = """samtools view -bS -"""

    child = subprocess.Popen(str(cline),
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=(sys.platform!="win32"),
                             close_fds = True)
    sout, serr = child.communicate(samdat)

    cline = """samtools sort - -o %s"""%(name)
    child2 = subprocess.Popen(str(cline),
                              stdin = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              stdout = subprocess.PIPE,
                              shell = (sys.platform!="win32"),
                              close_fds = True)
    bamdat, berr = child2.communicate(sout)

    with open(bamout, "wb") as o:
        o.write(bamdat)

    os.system("""samtools index %s > /dev/null 2> /dev/null"""%(bamout))
    bamidxdat = open(bamidx, "rb").read()

    if store:
        return samdat, bamdat, bamidxdat, bamout, bamidx
    else:
        return None, None, None, bamout, bamidx



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
    samdat, bamdat, bamidx, bamname, bamidxname = buildBAM(prefix, samname, False)
    

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
