#!/usr/bin/env python2.7

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


#http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks-in-python
def chunks(l, n):
    """ 
    Yield successive n-sized chunks from l. 
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def concatAlignment(con, specnumber, outfile):
    tally = [ row[0] for row in con.execute("""SELECT fileID, count(*) as 'size' FROM groups GROUP BY fileID HAVING size = ?;""", (specnumber,) ) ]
    
    #TODO: Need a clean way to gather what species identifiers exist.
    species = set()
    dat = []
    payload = defaultdict(list)
    query = """SELECT A.IDs, A.SEQS FROM ( %(csr_as_seqeuences)s ) AS A  where A.fileID IN (""" % dict(csr_as_seqeuences = QUERY_TRIMMED_CSR_AS_SEQS)
    for ctally in chunks(tally, 100):
        for r in con.execute(query + ",".join("?"*len(ctally)) + """) ORDER BY A.fileID;""", tuple(ctally) ):        
            idents = [ s.split("_")[0] for s in r[0].split("\t") ]
            species.update( idents )
            dat.append( [ (i, s,)  for i, s in itertools.izip(idents, r[1].split("\t") ) ] )
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
    os.system("""seqret  -sequence %s.fasta -sformat1 fasta -outseq %s.msa -osformat2 clustal"""%(outfile, outfile) )


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
    
