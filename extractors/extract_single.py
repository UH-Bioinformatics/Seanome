#!/usr/bin/env python2.7                                                                                                                                                                         
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

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--id', type = str, help = "Either a fileID or filename" )
parser.add_argument('-d', '--database',  required = True, help = "Seanome sqlite db")

args = parser.parse_args()

con = sqlite3.connect(args.database, check_same_thread=False)
curs = con.cursor()

curs.execute("""SELECT C.name, A.bam, A.bamidx, B.sequence, A.sam                                                                                                                   
                FROM trimmed_inferSAM AS A
                JOIN trimmed_consensus AS B ON (A.fileID = B.fileID) JOIN files AS C ON(C.id = A.fileID)
                WHERE C.name = ? OR A.fileID = ?;""", (str(args.id), str(args.id)))

rows = [(r[0], r[3], bytearray(r[1]), bytearray(r[2]), r[4]) for r in curs ]
with open("%s.bam"%(rows[0][0]), "w") as o :
    o.write(rows[0][2])
with open("%s.bam.bai"%(rows[0][0]), "w") as o :
    o.write(rows[0][3])
with open("%s.con"%(rows[0][0]), "w") as o :
    o.write(">Consensus\n%s"%(rows[0][1]))
with open("%s.sam"%(rows[0][0]), "w") as o :
    o.write(rows[0][4])
#print rows
