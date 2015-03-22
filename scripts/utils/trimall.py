import os
import sys
import subprocess
import itertools
import sqlite3
import cStringIO as StringIO
from Bio import SeqIO

from threadpool import ProducerConsumer
from utils import removeFiles, CONSENSUS_NAME
from sqlitedb import QUERY_CSR_AS_SEQS


def producer(info):
    inputdata = "%s.fasta"%(str(info[0]))
    consensus = str(info[1])
    seqs = [">%s\n%s"%(str(i), str(s)) for i, s in zip(info[2].split("\t"), info[3].split("\t"))]
    with open(inputdata, "w") as o:
        print >> o, ">%(seqID)s\n%(seq)s"%(CONSENSUS_NAME, consensus)
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


def runTrimAL(args, con):
    curs = con.cursor()
    curs.execute("""SELECT A.fileID, B.sequence, A.IDs, A.SEQS 
                    FROM ( %(csr_as_seqs)s ) AS A
                    JOIN consensus AS B ON (A.fileID = B.fileID);"""%dict(csr_as_seqs = QUERY_CSR_AS_SEQS ) )
    worker = ProducerConsumer(args, args.threads, producer, consumer)
    worker.run(con, curs)
