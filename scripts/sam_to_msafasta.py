#!/usr/bin/env python2.7

import re
import sys
import os
import argparse
import sqlite3
from Bio import SeqIO
from Bio.Seq import Seq
import pysam 
from utils.cigar import expandCigar, compressCigar
from utils.utils import removeFiles, CONSENSUS_NAME


def expandSAM(args, curs):
    for row in curs:
        
        bamfile = os.path.join(args.work, "%s_%s.bam"%(args.prefix, row[0]))
        bamidxfile = os.path.join(args.work, "%s_%s.bam.bai"%(args.prefix, row[0]))
        with open(bamfile, "wb") as o:
            o.write(row[1])
        with open(bamidxfile, "wb") as o:
            o.write(row[2])
        sfile = pysam.Samfile(bamfile)
        print ">%s\n%s"%(CONSENSUS_NAME, row[3])
        for read in sfile:
            newseq = "-"*(read.pos)
            cigar = expandCigar(read.cigarstring)
            seq = list(read.query)
            idx = 0
            for c in cigar:
                if c == 'M':
                    newseq += seq[idx]
                    idx += 1
                else:
                    newseq += '-'
            print ">%s\n%s"%(read.qname, newseq)
        removeFiles([bamfile, bamidxfile])
        


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-w', '--work', type = str, default = ".", help = "Working directory path" )
    parser.add_argument('-i', '--id', type = str, required = True, help = "Either a fileID or filename" )
    parser.add_argument('-d', '--database',  required = True, help = "Seanome sqlite db")
    parser.add_argument('-t', '--trimmed',  required = False, action = "store_true", help = "Generate information from the trimmed results")
    parser.add_argument('-p', '--prefix',  required = True, default = "", help = "Prefix for temp files")
    
    args = parser.parse_args()

    con = sqlite3.connect(args.database, check_same_thread=False)
    curs = con.cursor()
    if args.trimmed:
        curs.execute("""SELECT C.name, A.bam, A.bamidx, B.sequence                                                                                                                   
                FROM trimmed_inferSAM AS A
                JOIN trimmed_consensus AS B ON (A.fileID = B.fileID) JOIN files AS C ON(C.id = A.fileID)
                WHERE C.name = ? OR A.fileID = ?;""", (str(args.id), str(args.id)))
    else:
        curs.execute("""SELECT C.name, A.bam, A.bamidx, B.sequence
                FROM inferSAM AS A
                JOIN consensus AS B ON (A.fileID = B.fileID) JOIN files AS C ON(C.id = A.fileID)
                WHERE C.name = ? OR A.fileID = ?;""", (str(args.id), str(args.id)))


    expandSAM(args, curs)

    con.commit()
    con.close()
