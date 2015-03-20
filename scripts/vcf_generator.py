#!/usr/bin/python 
import argparse
import os
import sys
import subprocess
from multiprocessing import Pool
import sqlite3
import cStringIO as StringIO

from utils.threadpool import ProducerConsumer
from utils.utils import removeFiles

def producer(info):
    fileidx = str(info[0])

    inputdata = "%s.con.fasta"%(fileidx)
    consensus = str(info[1])

    with open(inputdata, "w") as o:
        print >> o, ">Consensus\n%s"%(consensus)
    os.system("""samtools faidx %s"""%(inputdata))

    bamfile = "%s.bam"%(fileidx)

    for g in info[3].split("\t"):
        cline = """samtools view -bhr "%s" - > %s_%s"""%(g, g, bamfile)
        child = subprocess.Popen(str(cline),
                                 stdin=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=(sys.platform!="win32"),
                                 close_fds = True)
        child.communicate(info[2])
        with open("%s.ids"%(fileidx), "a") as o:
            print >> o, "%s_%s"%(g, bamfile)   

    cline = """varscan.sh %s %s %s 2> /dev/null"""%( "%s.ids"%(fileidx), bamfile, inputdata )
    child = subprocess.Popen(str(cline),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=(sys.platform!="win32"),
                             close_fds = True)
    dat, err = child.communicate()
    removeFiles([inputdata, "%s.fai"%(inputdata)])
    return fileidx, dat


def consumer(con, returndata):
    curs = con.cursor()
    for r in returndata:
        if not r:
            continue
        fileid = int(r[0])
        curs.execute("""INSERT INTO trimmed_vcf(fileID, vcf) VALUES(?, ?)""", (fileid, r[1], ) )
    con.commit()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--threads', type = int, default = 1, help = "Number of processing threads. (default: 1)" )
    parser.add_argument('-d', '--database',  required = True, help = "Seanome sqlite db")

    args = parser.parse_args()

    con = sqlite3.connect(args.database, check_same_thread=False)
    con.execute("""PRAGMA foreign_keys = ON;""")
    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_vcf(id INTEGER PRIMARY KEY, fileID INTEGER, vcf TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_vcf_fileid_idx ON trimmed_vcf(fileID ASC);""")
    curs = con.cursor()

    curs.execute("""SELECT A.fileID, A.bam, B.sequence, C.groupids
                    FROM trimmed_inferSAM AS A 
                    JOIN trimmed_consensus AS B ON (A.fileID = B.fileID) 
                    JOIN (SELECT  D.fileID as fid, group_concat(D.groupid, '\t') as groupids FROM groups AS D GROUP BY fileID) AS C  ON (A.fileID = C.fid);""")

    rows = ( (r[0], r[2], bytearray(r[1]), r[3]) for r in curs )

    worker = ProducerConsumer(args, args.threads, producer, consumer)
    worker.run(con, rows)
    con.commit()
    con.close()
