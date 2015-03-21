import os
import sys
import subprocess
import itertools
import sqlite3
import cStringIO as StringIO
from Bio import SeqIO

from threadpool import ProducerConsumer
from utils import removeFiles


def producer(info):
    inputdata = "%s.fasta"%(str(info[0]))
    consensus = str(info[1])
    seqs = [">%s\n%s"%(str(i), str(s)) for i, s in zip(info[2].split("\t"), info[3].split("\t"))]
    with open(inputdata, "w") as o:
        print >> o, ">Consensus\n%s"%(consensus)
        print >> o, "%s"%("\n".join(seqs))
    cline = """trimal -in %s  -fasta -gt 0.8 -st 0.001 -cons 60 -colnumbering"""%(inputdata)
    child = subprocess.Popen(str(cline),
                             stdout=subprocess.PIPE,
                             universal_newlines = True,
                             shell=(sys.platform!="win32"))
    sout, serr = child.communicate()
    removeFiles([inputdata])
    sout = filter(None, sout.splitlines()) # strip empty lines
    fasta = "\n".join(sout[:-1])
    log = sout[-1]
    return fasta, log, info[0]


def consumer(con, returndata):
    curs = con.cursor()
    for r in returndata:
        if not r:
            continue
        seqs = list(SeqIO.parse(StringIO.StringIO(r[0]), "fasta"))
        curs.executemany("""INSERT INTO trimmed_csr(fileID, seqID, sequence) VALUES(?, ?, ?)""", 
                         [(idx, s.id, str(s.seq)) for idx, s in itertools.izip( itertools.repeat(r[2]), seqs[1:]) ])
        curs.execute("""INSERT INTO trimmed_consensus(fileID, sequence) VALUES(?, ?)""", (r[2], str(seqs[0].seq)))
        curs.execute("""INSERT INTO trimmed_logs(fileID, positions) VALUES(?, ?)""", (r[2], r[1]))
    con.commit()




def runTrimAL(args)
    con = sqlite3.connect(args.database, check_same_thread=False)
    con.execute("""PRAGMA foreign_keys = ON;""")
    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_csr(id INTEGER PRIMARY KEY, fileID INTEGER, seqID TEXT, sequence TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_csr_fileid_idx ON trimmed_csr(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_consensus(id INTEGER PRIMARY KEY, fileID INTEGER, sequence TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_con_fileid_idx ON trimmed_consensus(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_logs(id INTEGER PRIMARY KEY, fileID INTEGER, positions TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_log_fileid_idx ON trimmed_logs(fileID ASC);""")

    curs = con.cursor()
    curs.execute("""SELECT A.fileID, B.sequence, A.IDs, A.SEQS 
                    FROM (
                          SELECT fileID, group_concat(seqID, '\t') AS IDs, group_concat(sequence, '\t') AS SEQS
                          FROM csr group by fileID) AS A 
                    JOIN consensus AS B ON (A.fileID = B.fileID);""")

    worker = ProducerConsumer(args, args.threads, producer, consumer)
    worker.run(con, curs)
    con.commit()
    con.close()
