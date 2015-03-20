#!/usr/bin/python

import os
import sys
import argparse
import logging
import pysam
import tempfile
import shutil
import stat
import itertools
import sqlite3
import subprocess
import re

from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO

from utils.inferSAM import SAM_BUILDER
from utils.sqlitedb import buildsqlitedb
from utils.threadpool import ProducerConsumer


hmmer_hit_re = re.compile(r'score\s*bias\s*Evalue\s*hmmfrom\s*hmm\s*to\s*alifrom')
hmmer_query_length = re.compile("Query:\s*ali\s*\[M=(\d+)")


#TODO
# - ADD POSSIBILITY TO WORK ON A SINGLE FILE, RATHER THAN A DIRECTORY  


#https://stackoverflow.com/questions/1889597/deleting-directory-in-python/1889686#1889686
# to handle the case of read only files in the remove tree for shutil.rmtree
def remove_readonly(func, path, excinfo):
    os.chmod(path, stat.S_IWRITE)
    func(path)

def runInstance(myInstance):
   dryRunInstance(myInstance)
   myInstance.run()


def dryRunInstance(myInstance):
   logging.warning(myInstance.dryRun())


def makeDirOrdie(dirPath):
   if not os.path.isdir(dirPath):
      os.makedirs(dirPath)
   else:
       # TODO; Uncomment after testing done
      #logging.error("Split fasta directory %s already exists " % hmmerOutputDir)
      #sys.exit()
      pass
   return dirPath



def find_seed_csr(args):
    logging.debug("Finding seed common shared regions in %s: \n Starting... "%(args.input))

    inFile = open(args.input, 'r')
    outputdb= args.database
    outFileNum=0
    con = buildsqlitedb(args, outputdb)
    filemap = dict()
    c = con.cursor()
    for line in inFile:
        data = line.rstrip().split()
        aliLength = min(len(data[11]),len(data[12]))
        simRatio = float(data[2])
        if aliLength >= args.min_csr_len and simRatio >= args.min_csr_sim:
            fname = "outFile_%s_%s"%(str(outFileNum), str(aliLength))
            if fname not in filemap:
                c.execute("""INSERT INTO files(name) values(?)""", (fname,) )
                c.execute("""SELECT id FROM files WHERE name = ?;""", (fname,))
                row = c.fetchone()
                filemap[fname] = row[0]
                fileID = filemap[fname]
                seqid1 = "%s_%s_%s%s"%(data[0], data[7], data[8] , "_R" if(data[5] == "-") else "" )
                seqid2 = "%s_%s_%s%s"%(data[1], data[9], data[10], "_R" if(data[6] == "-") else "")
            c.execute("""INSERT INTO csr(fileID, seqID, sequence, seed) values(?,?,?,?)""", (fileID, seqid1, data[11], 1) )
            c.execute("""INSERT INTO csr(fileID, seqID, sequence, seed) values(?,?,?,?)""", (fileID, seqid2, data[12], 1) )
            outFileNum += 1

    con.execute("""CREATE INDEX IF NOT EXISTS csr_fileid_idx ON csr(fileID ASC);""")
    con.execute("""CREATE INDEX IF NOT EXISTS csr_seed_idx ON csr(seed ASC);""")
    con.commit()
    con.close()
    inFile.close()
    logging.debug("Finished finding seed common shared regions in %s"%(args.input))




def getSubSeqFromFasta(seq, start, end, reverse=False):
    ident = "%s_%s_%s%s" % (seq.id, start + 1, end, "_R" if(reverse) else "" )
    if reverse:
        seqo = seq[start:end].reverse_complement().seq
    else:
        seqo = seq[start:end].seq
    return ident, seqo


def find_csrProducer(payload):
    info, genome, genomeseq, min_csr_len, min_csr_sim = payload
    seqs = "\n".join([">%s\n%s"%(str(i), str(s)) for i, s in zip(info[1].split("\t"), info[2].split("\t"))])
    cline = """nhmmer --cpu 1 --qformat fasta - %s"""%(genome)
    child = subprocess.Popen(str(cline),
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             universal_newlines=True,
                             shell=(sys.platform!="win32"),
                             close_fds = True)
    hitfile, hiterr = child.communicate(seqs)
    for line in hitfile.splitlines():
        line = line.strip()
        if hmmer_query_length.search(line):
            length = hmmer_query_length.search(line).groups()[0]
        elif hmmer_hit_re.search(line):
            line = hitfile.next()
            line = hitfile.next()
            data = line.strip().split()
            qStart, qEnd, hStart, hEnd, acc = [int(data[i]) for i in [4, 5, 10, 11]] + [ data[14] ]
            seqLen = 0
            if  hEnd > hStart:
                strand="+"
                seqLen = hEnd - hStart + 1
            else:
                seqLen = hStart - hEnd + 1
                strand="-"

            if int(seqLen) > min_csr_len and float(acc) > min_csr_sim:
                if strand=="+":
                    hitid, hitseq = getSubSeqFromFasta(genomeseq, hStart - 1, hEnd) # nhmmer alignment is 1 based 
                else:
                    hitid, hitseq = getSubSeqFromFasta(genomeseq, hEnd - 1, hStart, True) # nhmmer alignment is 1 based 
            else:
                print "Parameters not met for alignment %s " % info[0]
                continuehitid
            # need to run muscle on the seqs + hitseq
            seqs = "\n".join([">%s\n%s"%(str(i), str(s)) for i, s in zip(info[1].split("\t"), info[2].replace("-","").split("\t"))])            
            seqs = ">%s\n%s\n%s"%(hitid, hitseq, seqs)
            cline = MuscleCommandline(quiet = True)
            childM = subprocess.Popen(str(cline),
                                      stdin=subprocess.PIPE,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE,
                                      universal_newlines=True,
                                      shell=(sys.platform!="win32"),
                                      close_fds = True)
            mucout, mucerr = childM.communicate(seqs)
            alnA = AlignIO.read(StringIO.StringIO(mucout), "fasta")
            return info[0], alnA, hitid
    return None

    
def find_csrConsumer(con, returndat):
   cur2 = con.cursor()
   for dat in returndat:
       if not dat:
           continue
       print dat
       for q in dat[1]:
           print str(q)
           if q.id == dat[2]:
               cur2.execute("""INSERT INTO csr(fileID, seqID, sequence) VALUES(?, ?, ?)""", (dat[0], q.id, str(q.seq), ) )
           else:
               cur2.execute("""UPDATE csr SET sequence = ?  WHERE seqID = ?;""", (str(q.seq), q.id,) )
   con.commit()


def find_csr(args):
   con = sqlite3.connect(args.database, check_same_thread=False)
   cur = con.cursor()
   cur.execute("""SELECT fileID, group_concat(seqID, '\t') AS IDs, group_concat(sequence, '\t') AS SEQS 
                   FROM  csr 
                   GROUP BY fileID;""")
   genomeseq = SeqIO.read(args.genome, "fasta")

   tmpDir = tempfile.mkdtemp(prefix="seanome_", dir=".")
   nhmmer_dir = os.path.join(tmpDir, "nhmmer_dir")
   os.mkdir(nhmmer_dir)
   worker = ProducerConsumer(args, args.threads, find_csrProducer, find_csrConsumer)
   worker.run(con, ( (c, args.genome, genomeseq, args.min_csr_len, args.min_csr_sim) for c in cur) )
   con.commit()
   con.close()
   # DLS cleanup after mkdtemp.  It is the Users responsibility!
   try:
      shutil.rmtree(tmpDir, onerror = remove_readonly)
   except:
      pass


def inferSAMProducer(info):
    #print info
    name = str(info[0])
    fastseqs = info[1] + info[2]
    samdir = info[3]
    sb = SAM_BUILDER(name, fastseqs, samdir)
    return sb.run()


def inferSAMConsumer(con, returndat):
    curs = con.cursor()
    for dat in ret:
        if not dat:
            continue
        sam = dat[0]
        bam = dat[1]
        bamidx = dat[2]
        ident = dat[3]
        groups = dat[4]
        newcon = dat[5]
        curs.execute("""INSERT INTO inferSAM(fileID, sam, bam, bamidx) VALUES(?,?,?,?);""", (ident, sam, sqlite3.Binary(bam), sqlite3.Binary(bamidx),))
        curs.execute("""UPDATE consensus SET sequence = ? WHERE fileID = ?;""", (newcon, ident, ) )
        curs.executemany("""INSERT INTO groups(fileID, groupid) VALUES(?,?);""", ( (ident, g,) for g in groups )  )
    con.commit()

def inferSAM(args):
    con = sqlite3.connect(args.database, check_same_thread=False)
    con.execute("""PRAGMA foreign_keys = ON;""")
    con.execute("""CREATE TABLE IF NOT EXISTS inferSAM(id INTEGER PRIMARY KEY, fileID INTEGER, sam TEXT, bam BLOB, bamidx BLOB, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS infersam_fileid_idx ON inferSAM(fileID ASC);""")
    con.execute("""CREATE TABLE IF NOT EXISTS groups(id INTEGER PRIMARY KEY, fileID INTEGER, groupid TEXT, FOREIGN KEY(fileID) REFERENCES files(id) );""")
    con.execute("""CREATE INDEX IF NOT EXISTS groups_fileid_idx ON groups(fileID ASC);""")

    con.commit()

    curs = con.cursor()
    curs.execute("""SELECT A.fileID, B.sequence, A.IDs, A.SEQS 
                    FROM (
                          SELECT fileID, group_concat(seqID, '\t') AS IDs, group_concat(sequence, '\t') AS SEQS 
                          FROM csr group by fileID) AS A 
                    JOIN consensus AS B ON (A.fileID = B.fileID);""")
    rows = (  (r[0], ">Consensus\n%s\n"%(r[1]), "\n".join([">%s\n%s"%(i,s,) for i, s in itertools.izip(r[2].split("\t") , r[3].split("\t") ) ] ) , sd,) for r, sd in itertools.izip(curs, itertools.repeat(str(args.split_sam_dir)) )    ) 

    worker = ProducerConsumer(args, args.threads, inferSAMProducer, inferSAMConsumer)
    worker.run(con, rows)
    con.commit()
    con.close()


def main(argv):
    parser = argparse.ArgumentParser(description="Seanome description", epilog="Seanome long text description")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s '+version)
    parser.add_argument('-t', '--threads', type=int, default = 1)
    subparsers = parser.add_subparsers(dest='action', help='Available commands')

    # seed_csr
    parser_mask = subparsers.add_parser('seed_csr')
    parser_mask.add_argument('-i', '--input',  required=True, help="Lastz like input file")
    parser_mask.add_argument('-d', '--database', required=True, help="Seanome sqlite database")
    parser_mask.add_argument('-l', '--min_csr_len', type=int, default=150, help=" Minimum common shared region length")
    parser_mask.add_argument('-s', '--min_csr_sim', type=float, default=0.88, help=" Minimum common shared region similarity")
    parser_mask.set_defaults(func = find_seed_csr)

    # csr
    parser_mask = subparsers.add_parser('find_csr')
    parser_mask.add_argument('-d', '--database', required=True, help="Seanome sqlite database")
    parser_mask.add_argument('-g', '--genome', required=True, help="Reference genome") # pseudo genome 
    parser_mask.add_argument('-l', '--min_csr_len', type=int, default=150, help=" Minimum common shared region length")
    parser_mask.add_argument('-s', '--min_csr_sim', type=float, default=0.88, help=" Minimum common shared region similarity")
    parser_mask.set_defaults(func = find_csr)

    # Build transitive SAM from MSA
    parser_mask = subparsers.add_parser('inferSAM')
    parser_mask.add_argument('-d', '--database', required=True, help="Seanome sqlite database")
    parser_mask.add_argument('-s', '--split_sam_dir', required = True, help = 'Split SAM Files directory')
    parser_mask.set_defaults(func=inferSAM)

    # consensus
    #parser_mask = subparsers.add_parser('consensus')
    #parser_mask.add_argument('-a', '--alis_dir',  required=True, help="Inut alignments directory")
    #parser_mask.add_argument('-c', '--cons_output', required=True, type=makeDirOrdie, help="Output Consensus Directory")
    #parser_mask.set_defaults(func=consensus)

    # # Mask Genome
    # parser_mask = subparsers.add_parser('mask')
    # parser_mask.add_argument('-i', '--input',  required=True, help=" Input file to maks")
    # parser_mask.add_argument('-o', '--output', required=True, help=" Masked output file")
    # parser_mask.add_argument('-l', '--min_length', default=80, help=" Minimum alignment length")
    # parser_mask.add_argument('-s', '--min_similarity', default=0.86, help=" Minimum alignment similarity")
    # parser_mask.set_defaults(func=maskGenome)

    # Parse arguments
    args = parser.parse_args()
    #pool = Pool(processes=args.threads)
    #logging.debug("Initial ARGS are:")
    #logging.debug(args)
    args.func(args)


if __name__ == "__main__":
   # test
   version = "alpha 0.01"
   FORMAT = "%(asctime)-15s  %(message)s"
   logging.basicConfig(format=FORMAT, level=logging.DEBUG)
   main(sys.argv)
