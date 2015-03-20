#!/usr/bin/python

import sys
import sqlite3
import argparse
import itertools
import cStringIO as StringIO
from multiprocessing import Pool
from Bio import AlignIO
from Bio.Align import AlignInfo


def processor(info):
    fid = info[0]
    seqs = "\n".join([">%s\n%s"%(str(i), str(s)) for i, s in zip(info[1].split("\t"), info[2].split("\t"))])
    i = AlignIO.read(StringIO.StringIO(seqs), 'fasta')        
    si = AlignInfo.SummaryInfo(i)
    mySeq = si.dumb_consensus(threshold = 0.1, ambiguous = "N")
    return info[0], str(mySeq)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--threads', type = int, default = 1, help = "Number of processing threads. (default: 1)" )
    parser.add_argument('-d', '--database',  required = True, help = "Seanome sqlite db")
    args = parser.parse_args()

    pool = Pool(processes = args.threads)
    con = sqlite3.connect(args.database, check_same_thread = False)
    con.execute("""PRAGMA foreign_keys = ON;""")
    con.execute("""CREATE TABLE IF NOT EXISTS consensus(id INTEGER PRIMARY KEY, fileID INTEGER, sequence TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS con_fileid_idx ON consensus(fileID ASC);""")
    con.commit()
    curs = con.cursor()
    curs2 = con.cursor()
    
    curs.execute("""SELECT fileID, group_concat(seqID, '\t') AS IDs, group_concat(sequence, '\t') AS SEQS FROM csr group by fileID;""")
    ret = pool.imap_unordered(processor, curs)
    for dat in ret:
        if not dat:
            continue
        curs2.execute("""INSERT INTO consensus(fileID, sequence) VALUES(?, ?)""", (dat[0], dat[1],) )
    con.commit()
    con.close()
