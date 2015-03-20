#!/usr/bin/python 
import argparse
import os
import sys
import subprocess
import itertools
from multiprocessing import Pool

import sqlite3
import cStringIO as StringIO


from Bio import SeqIO
from Bio import AlignIO

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio.Align.Applications import MuscleCommandline


def trimalprocess(info):
    """
    A thread pool worker function.  
    - It will run muscle using stdin and stdout (no intermediate files required!
    - Build a consensus from the alignment Muscle finds
    - Drop ambigious clusters if requested
    - Build a new cigar based on the MSA alignment and the consensus that was created.
    - works on a single MSA and returns 1 consensus, with multiple cigards (1 for each sequence)
    """
#A.fileID, B.sequence, A.IDs, A.SEQS
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
    try:
        os.remove(inputdata)
    except:
        pass
    sout = filter(None, sout.splitlines())
    fasta = "\n".join(sout[:-1])
    log = sout[-1]
    return fasta, log, info[0]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--threads', type = int, default = 1, help = "Number of processing threads. (default: 1)" )
    parser.add_argument('-d', '--database',  required = True, help = "Seanome sqlite db")

    args = parser.parse_args()

    pool = Pool(processes = args.threads)

    con = sqlite3.connect(args.database, check_same_thread=False)
    con.execute("""PRAGMA foreign_keys = ON;""")
    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_csr(id INTEGER PRIMARY KEY, fileID INTEGER, seqID TEXT, sequence TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_csr_fileid_idx ON trimmed_csr(fileID ASC);""")
    #con.execute("""CREATE INDEX IF NOT EXISTS trimmed_csr_seed_idx ON trimmed_csr(seed ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_consensus(id INTEGER PRIMARY KEY, fileID INTEGER, sequence TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_con_fileid_idx ON trimmed_consensus(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_logs(id INTEGER PRIMARY KEY, fileID INTEGER, positions TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_log_fileid_idx ON trimmed_logs(fileID ASC);""")


    curs = con.cursor()


    #curs.execute("""SELECT C.name, A.fileID, B.sequence, A.IDs, A.SEQS FROM (SELECT fileID, group_concat(seqID, '\t') AS IDs, group_concat(sequence, '\t') AS SEQS FROM sharedparts group by fileID) AS A JOIN consensus AS B ON (A.fileID = B.fileID) JOIN files AS C ON (C.id = A.fileID);""")
    curs.execute("""SELECT A.fileID, B.sequence, A.IDs, A.SEQS 
                    FROM (
                          SELECT fileID, group_concat(seqID, '\t') AS IDs, group_concat(sequence, '\t') AS SEQS 
                          FROM csr group by fileID) AS A 
                    JOIN consensus AS B ON (A.fileID = B.fileID);""")
    #rows = ( (r[0], r[1], r[2], r[3]) for r in curs )
    ret = pool.imap_unordered( trimalprocess, curs )

    for r in ret:
        if not r:
            continue
        seqs = list(SeqIO.parse(StringIO.StringIO(r[0]), "fasta"))
        con.executemany("""INSERT INTO trimmed_csr(fileID, seqID, sequence) VALUES(?, ?, ?)""", 
                      [(idx, s.id, str(s.seq),) for idx, s in itertools.izip( itertools.repeat(r[-1]), seqs[1:]) ])
        con.execute("""INSERT INTO trimmed_consensus(fileID, sequence) VALUES(?, ?)""", (r[-1], str(seqs[0].seq)))
        con.execute("""INSERT INTO trimmed_logs(fileID, positions) VALUES(?, ?)""", (r[-1], r[1]))
        
    pool.close()
    con.commit()
    con.close()
