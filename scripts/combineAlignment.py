#!/usr/bin/python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sqlite3
import cStringIO as StringIO
import argparse
import os
from utils.sqlitedb import QUERY_CSR_AS_SEQS, QUERY_TRIMMED_CSR_AS_SEQS
from collections import defaultdict
import itertools 

def concatAlignment(con, specnumber, outfile):
    tally = [ row[0] for row in con.execute("""SELECT A.fileID FROM (SELECT fileID, count(*) as 'size' FROM trimmed_csr GROUP BY fileID ) AS A WHERE size = ? ;""", (specnumber,) ) ]
    
    #TODO: Need a clean way to gather what species identifiers exist.

    species = set()
    dat = []
    payload = defaultdict(list)

    for r in con.execute("""SELECT A.IDs, A.SEQS FROM ( %(csr_as_seqeuences)s ) AS A  where A.fileID IN ("""%dict(csr_as_seqeuences = QUERY_TRIMMED_CSR_AS_SEQS) + ",".join("?"*len(tally)) + """) ORDER BY A.fileID;""", tuple(tally) ):        
        species.update( ( s.split("_")[0] for s in r[0].split("\t") ) )
        dat.append( [ (i.split("_")[0],s,)  for i, s in itertools.izip(r[0].split("\t") , r[1].split("\t") ) ] )
    species = list(species)

    for d in dat:
        notused = dict.fromkeys(species)
        slen = 0
        for ident, seq in d:
            payload[ident].append(seq)
            slen = len(seq)
            notused.pop(ident, None)
        for k in notused.iterkeys():
            payload[k].append("-"*slen)
    seqs = []
    with open("%s.fasta"%(outfile), "w") as o:
        for k, v in payload.iteritems():
            s = "".join(v)
            print >> o, ">%s\n%s"%(k, s)
            seqs.append(SeqRecord(Seq(s), id = k) )
    SeqIO.write(seqs, "%s.msa"%(outfile), "clustal")




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--count', type = int, required = True, help = "Number of Species involved" )
    parser.add_argument('-d', '--database',  required = True, help = "Seanome sqlite db")
    parser.add_argument('-o', '--output',  required = True, help = "name of the output file")
    args = parser.parse_args()

    con = sqlite3.connect(args.database, check_same_thread=False)
    concatAlignment(con, args.count, args.output)
    con.commit()
    con.close()
    
