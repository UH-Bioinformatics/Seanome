#!/usr/bin/python

from Bio import SeqIO
import pysam 
import sys
import re
import sqlite3
import cStringIO as StringIO
import argparse
import os
import sys
import subprocess
import itertools
from multiprocessing import Pool


cigar_re = re.compile(r'([0-9]+)([M=XID])')
def expandCigar( cigar):
    cigstr = ""
    for m in cigar_re.finditer(cigar):
        cigstr +=  m.group(2) * int(m.group(1))
    return cigstr


def compressCigar(cigar):
    cigar = cigar.rstrip("D")
    c = None
    cnt = 0
    ciglst = []
    for b in cigar:
        if b != c:
            if c!= None:
                ciglst.append("%s%s"%(cnt, c))
            cnt = 1
            c = b
        else:
            cnt += 1
    if c!= None:
        ciglst.append("%s%s"%(cnt, c))
    return "".join(ciglst)


def processor(info):
    fileidx = str(info[0])   
    consensus = ">Consensus\n%s\n"%(str(info[2]))    
    seqs = [">%s\n%s"%(str(i), str(s)) for i, s in zip(info[3].split("\t"), info[4].split("\t"))]
    
    fastafile = StringIO.StringIO(consensus + "\n".join(seqs) )
    
    bamfile = "%s.bam"%(fileidx)
    bamidxfile = "%s.bam.bai"%(fileidx)
    with open(bamfile, "wb") as o:
        o.write(info[5])
    with open(bamidxfile, "wb") as o:
        o.write(info[6])

    refs = [ (r.id, len(r.seq) ) for r in  SeqIO.parse(fastafile, "fasta")]

    sfile = pysam.Samfile(bamfile)

    hdr = sfile.header.copy()
    hdr['SQ'] = [ {'LN': refs[0][1], 'SN': refs[0][0] }]
    
    total = 0
    indices = eval("[" + info[1] + "]")
    maxidx = max(indices) + 1
    shift = [0]*(maxidx)
    for i in xrange(len(shift)):
        shift[i] = total
        if i not in indices:
            total += 1
    outfile = pysam.Samfile("%s.trim.sam"%(fileidx), "wh", header = hdr)
    for read in sfile.fetch():
        newseq = ""
        newcig = ""
        newqual = ""
        cigar = expandCigar(read.cigarstring)
        back = 0
        pos = read.pos
        for idx, c in enumerate(cigar):
            # stop if we exceed the reference length
            if pos >= maxidx:
                break
            pos += 1
            # only process bases that we havent removed in cleanup
            # also, count any D's that we missed to shift us back in the query and qual
            if (idx + read.pos) not in indices:
                # an idx not in indices indicates that trimal has removed this particular position
                if cigar[idx].upper() == 'D':
                    back += 1
                continue

            # if the idx exists in the logs of trimAL, the base should be kept.
            newcig += cigar[idx]
            #D means we dont have a base to add.
            if cigar[idx].upper() != 'D': 
                newseq += read.query[idx - back]
                if read.qqual:
                    newqual += read.qqual[idx - back]
            else: 
                back += 1

        #if a sequence has nothing to add, why add it ?!?
        if newseq:
            read.pos = max( (read.pos - shift[read.pos]), 0)
            read.seq = newseq
            read.cigarstring = compressCigar(newcig)
            if not newqual:
                read.qual = "I"* len(newseq)
            else:
                read.qual = newqual
            outfile.write(read)
    outfile.close()
    
    samout = "%s.trim.sam"%(fileidx)
    bamout = "%s.trim.bam"%(fileidx)
    
    try:
        os.remove(bamfile)
    except:
        pass
    try:
        os.remove(bamidxfile)
    except:
        pass
    
    samdat = open(samout).read()
    try:
        os.remove(samout)
    except:
        pass

    cline = """samtools view -bS -"""
    child = subprocess.Popen(str(cline),
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=(sys.platform!="win32"),
                             close_fds = True)
    sout, serr = child.communicate(samdat)
    
    cline = """samtools sort - -o %s.trim"""%(fileidx)
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
    bamidxdat = open("%s.bai"%(bamout), "rb").read()
    try:
        os.remove(bamout)
    except:
        pass    
    try:
        os.remove("%s.trim.bam.bai"%(fileidx))
    except:
        pass    

    return fileidx, samdat, bamdat, bamidxdat
    #return fileidx, samout, bamout, "%s.trim.bam.bai"%(fileidx)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--threads', type = int, default = 1, help = "Number of processing threads. (default: 1)" )
    parser.add_argument('-d', '--database',  required = True, help = "Seanome sqlite db")
    args = parser.parse_args()

    pool = Pool(processes = args.threads)


    con = sqlite3.connect(args.database, check_same_thread=False)
    con.execute("""PRAGMA foreign_keys = ON;""")
    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_inferSAM(id INTEGER PRIMARY KEY, fileID INTEGER, sam TEXT, bam BLOB, bamidx BLOB, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_infersam_fileid_idx ON trimmed_inferSAM(fileID ASC);""")
    con.commit()

    curs = con.cursor()
    
    curs.execute("""SELECT A.fileID, A.positions, B.sequence, C.SEQS, C.IDs, D.bam, D.bamidx 
                    FROM trimmed_logs AS A 
                         JOIN trimmed_consensus AS B ON (A.fileID = B.fileID) 
                         JOIN (SELECT fileID, group_concat(seqID, '\t') AS IDs, group_concat(sequence, '\t') AS SEQS 
                               FROM trimmed_csr GROUP BY fileID) AS C ON (C.fileID = A.fileID) 
                         JOIN inferSAM AS D ON (D.fileID = A.fileID)""")  

    rows = ( (r[0], r[1], r[2], r[3], r[4], bytearray(r[5]), bytearray(r[6]),) for r in curs)   
    ret = pool.imap_unordered(processor, rows)
    curs2 = con.cursor() 
    for dat in ret:
        if not dat:
            continue
        sam = dat[1]
        bam = dat[2]
        bamidx = dat[3]
        ident = dat[0]
        curs2.execute("""INSERT INTO trimmed_inferSAM(fileID, sam, bam, bamidx) VALUES(?,?,?,?);""", (ident, sam, sqlite3.Binary(bam), sqlite3.Binary(bamidx),))
    con.commit()
    con.close()
    
