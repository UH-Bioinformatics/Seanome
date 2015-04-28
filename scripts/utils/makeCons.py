import sys
import sqlite3
import itertools
import cStringIO as StringIO

from Bio import AlignIO
from Bio.Align import AlignInfo

from sqlitedb import QUERY_CSR_AS_SEQS
from threadpool import ProducerConsumer


def producer(info):
    fid = info[0]
    seqs = "\n".join([">%s\n%s"%(str(i), str(s)) for i, s in zip(info[1].split("\t"), info[2].split("\t"))])
    i = AlignIO.read(StringIO.StringIO(seqs), 'fasta')        
    si = AlignInfo.SummaryInfo(i)
    mySeq = si.dumb_consensus(threshold = 0.0001, ambiguous = "N")
    return info[0], str(mySeq)


def consumer(con, returndata):
    curs = con.cursor()   
    for dat in returndata:
        if not dat:
            continue
        curs.execute("""INSERT INTO consensus(fileID, sequence) VALUES(?, ?)""", (dat[0], dat[1],) )
    con.commit()


def generateConsensusSequences(args, con):
    curs = con.cursor()    
    curs.execute("%s;"%(QUERY_CSR_AS_SEQS))
    worker = ProducerConsumer(args, args.threads, producer, consumer)
    worker.run(con, curs)
