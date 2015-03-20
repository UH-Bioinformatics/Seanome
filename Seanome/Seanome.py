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

from multiprocessing import Pool
from classes.ProgramRunner import *
from classes.NHMMER_TO_ALI import *
from classes.SAM_BUILDER import *


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


# def maskFile(args,samFile):
#    '''
#    Used in maskGenome
#    '''
#    samfile = pysam.Samfile(samFile)
#    alis = samfile.fetch()
#    mySeq =  SeqIO.read(args.input, 'fasta')
#    mseq = mySeq.seq.tomutable()
#    for ali in alis:
#       tags =  dict(ali.tags)
#         # if alignment hits with at least the minimum similarity, then mask the hit region                                                                                          
#       if tags['NM'] and (1 - float(tags['NM'])/int(ali.qlen) >= args.min_similarity):
#          # TODO: collect the intervals and do it at once since intervals overalap
#          maskRange = (min(ali.positions), max(ali.positions))
#          mseq[maskRange[0]: maskRange[1]+1] = "N" * (maskRange[1] - maskRange[0] + 1)
#          mySeq = SeqRecord(mseq, id=mySeq.id)
#    SeqIO.write(mySeq, args.output, 'fasta')


# def maskGenome(args, pool):
#    #print args.input
#    name =  os.path.splitext(os.path.basename(args.input))[0]    
#    # Can eventually be parallelized using the pool of threads
#    logging.debug("Starting the masking of input %s"%(args.input))
#    # run 
#    logging.debug("Computing frequencies for finding repeats" )
#    prog = ProgramRunner("build_lmer_table", [args.input, "%s.freqs"%(name) ] )

#    prog.run()
#    logging.debug("Finding repeats in the %s"%(args.input))
#    prog = ProgramRunner("repeatScout", [args.input, "%s_repeats.fa"%(name), name+".freqs" ] )

#    prog.run()
#    # Filter the repeats that do not pass the requirements
#    # for now, this only consists in dropping short reads < N
#    logging.debug("Filetering repeats")
#    repeats = SeqIO.parse("%s_repeats.fa"%(name), 'fasta')
#    filterReads=[]
#    for read in repeats:
#       if len(read.seq) >= args.min_length:
#          filterReads.append(read)
#    SeqIO.write(filterReads, "%s_repeats_filtered.fa"%(name), 'fasta')
#    if len(filterReads) >= 1:
#       prog = ProgramRunner("bowtie-build", [args.input, name])
#       prog.run()
#       prog = ProgramRunner("bowtie-align", [name, "%s_repeats_filtered.fa"%(name), "%s_repeats_filtered.sam"%(name)])
#       prog.run()
#       maskFile(args, "%s_repeats_filtered.sam"%(name))
#    else:
#       logging.debug("No repeats found is %s"%(args.input))
#    logging.debug("Done masking input %s"%(args.input))


def buildsqlitedb(args, dbname):
    con = sqlite3.connect(dbname)
    con.execute("""PRAGMA foreign_keys = ON;""")

    con.execute("""CREATE TABLE IF NOT EXISTS files(id INTEGER PRIMARY KEY, name TEXT);""")
    con.execute("""CREATE INDEX IF NOT EXISTS file_name_idx ON files(name ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS csr(id INTEGER PRIMARY KEY, fileID INTEGER, seqID TEXT, sequence TEXT, seed BOOLEAN DEFAULT 0, FOREIGN KEY(fileID) REFERENCES files(id));""")

    con.execute("""CREATE TABLE IF NOT EXISTS groups(id INTEGER PRIMARY KEY, fileID INTEGER, groupid TEXT, FOREIGN KEY(fileID) REFERENCES files(id) );""")
    con.execute("""CREATE INDEX IF NOT EXISTS groups_fileid_idx ON groups(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS consensus(id INTEGER PRIMARY KEY, fileID INTEGER, sequence TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS con_fileid_idx ON consensus(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_csr(id INTEGER PRIMARY KEY, fileID INTEGER, seqID TEXT, sequence TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_csr_fileid_idx ON trimmed_csr(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_consensus(id INTEGER PRIMARY KEY, fileID INTEGER, sequence TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_con_fileid_idx ON trimmed_consensus(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_logs(id INTEGER PRIMARY KEY, fileID INTEGER, positions TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_log_fileid_idx ON trimmed_logs(fileID ASC);""")
    
    con.execute("""CREATE TABLE IF NOT EXISTS inferSAM(id INTEGER PRIMARY KEY, fileID INTEGER, sam TEXT, bam BLOB, bamidx BLOB, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS infersam_fileid_idx ON inferSAM(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_inferSAM(id INTEGER PRIMARY KEY, fileID INTEGER, sam TEXT, bam BLOB, bamidx BLOB, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_infersam_fileid_idx ON trimmed_inferSAM(fileID ASC);""")

    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_vcf(id INTEGER PRIMARY KEY, fileID INTEGER, vcf TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_vcf_fileid_idx ON trimmed_vcf(fileID ASC);""")
 
    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_modvcf(id INTEGER PRIMARY KEY, fileID INTEGER, vcf TEXT, json TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_modvcf_fileid_idx ON trimmed_modvcf(fileID ASC);""")

    con.commit()
    return con

    
def generateAlis(args, lastzOutFile):
    ''' 
    Used in find_seed_csr 
    '''
    inFile = open(lastzOutFile, 'r')
    outputdb= args.database
    outFileNum=0
    con = buildsqlitedb(args, outputdb)
    filemap = dict()
    c = con.cursor()
    for line in inFile:
        line = line.rstrip()
        data = line.split()
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
            strand1 = ""
            strand2 = ""
            if data[5] == "-":
                strand1 = "_R"
            if data[6]  == "-":
                strand2 = "_R"
            c.execute("""INSERT INTO csr(fileID, seqID, sequence, seed) values(?,?,?,?)""", (fileID, "%s_%s_%s%s"%(data[0], data[7], data[8], strand1 ), data[11], 1) )
            c.execute("""INSERT INTO csr(fileID, seqID, sequence, seed) values(?,?,?,?)""", (fileID, "%s_%s_%s%s"%(data[1], data[9], data[10], strand2), data[12], 1) )
            outFileNum += 1

    con.execute("""CREATE INDEX IF NOT EXISTS csr_fileid_idx ON csr(fileID ASC);""")
    con.execute("""CREATE INDEX IF NOT EXISTS csr_seed_idx ON csr(seed ASC);""")
    con.commit()
    con.close()
    inFile.close()
    

def find_seed_csr(args, pool):
    # this process is not multi-threaded.. This might be somethign we want to do at a later stage?
    logging.debug("Finding seed common shared regions in %s: \n Starting... "%(args.input))
    generateAlis(args, args.input)
    logging.debug("Finished finding seed common shared regions in %s"%(args.input))


def getSubSeqFromFasta(seq, start, end, reverse=False):
    if reverse:
        ident = "%s_%s_%s_R" % (seq.id, start + 1, end)
        seqo = seq[start:end].reverse_complement().seq
    else:
        ident = "%s_%s_%s" % (seq.id, start + 1, end)
        seqo = seq[start:end].seq
    return ident, seqo


def process_find_csr(payload):
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
                continue
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

    
def find_csr(args, pool):
   # TODO
   # validate the params
   # throw error if program returns error

   tmpDir = tempfile.mkdtemp(prefix="seanome_", dir=".")
   nhmmer_dir = os.path.join(tmpDir, "nhmmer_dir")
   os.mkdir(nhmmer_dir)

   con = sqlite3.connect(args.database, check_same_thread=False)
   cur = con.cursor()
   cur.execute("""SELECT 
                         fileID, group_concat(seqID, '\t') AS IDs, group_concat(sequence, '\t') AS SEQS 
                   FROM 
                         csr 
                   GROUP BY
                         fileID;""")
   genomeseq = SeqIO.read(args.genome, "fasta")
   ret = pool.imap_unordered(process_find_csr, ( (c, args.genome, genomeseq, args.min_csr_len, args.min_csr_sim) for c in cur) )
   cur2 = con.cursor()
   for dat in ret:
       if not dat:
           continue
       print dat
       for q in dat[1]:
           print str(q)
           if q.id == dat[2]:
               cur2.execute("""INSERT INTO csr(fileID, seqID, sequence) VALUES(?, ?, ?)""", (dat[0], q.id, str(q.seq),) )
           else:
               cur2.execute("""UPDATE csr SET sequence = ?  WHERE seqID = ?;""", (str(q.seq), q.id,) )
   con.commit()
   con.close()
   # DLS cleanup after mkdtemp.  It is the Users responsibility!
   try:
      shutil.rmtree(tmpDir, onerror = remove_readonly)
   except:
      pass


def buildSam(info):
    #print info
    name = str(info[0])
    fastseqs = info[1] + info[2]
    samdir = info[3]
    sb = SAM_BUILDER(name, fastseqs, samdir)
    return sb.run()


def inferSAM(args, pool):
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
    
    ret = pool.imap_unordered( buildSam, rows)
    curs2 = con.cursor()
    for dat in ret:
        if not dat:
            continue
        sam = dat[0]
        bam = dat[1]
        bamidx = dat[2]
        ident = dat[3]
        groups = dat[4]
        newcon = dat[5]
        curs2.execute("""INSERT INTO inferSAM(fileID, sam, bam, bamidx) VALUES(?,?,?,?);""", (ident, sam, sqlite3.Binary(bam), sqlite3.Binary(bamidx),))
        curs2.execute("""UPDATE consensus SET sequence = ? WHERE fileID = ?;""", (newcon, ident, ) )
        curs2.executemany("""INSERT INTO groups(fileID, groupid) VALUES(?,?);""", ( (ident, g,) for g in groups )  )
    con.commit()
    con.close()
   

#def consensus(args, pool):
#   aliFiles = os.listdir(args.alis_dir)
#   pool.map(runInstance, [ProgramRunner("addConsensus", [os.path.join(args.alis_dir, x), os.path.join(args.cons_output, x)] ) for x in aliFiles])


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
    pool = Pool(processes=args.threads)
    #logging.debug("Initial ARGS are:")
    #logging.debug(args)
    args.func(args, pool)


if __name__ == "__main__":
   # test
   version = "alpha 0.01"
   FORMAT = "%(asctime)-15s  %(message)s"
   logging.basicConfig(format=FORMAT, level=logging.DEBUG)
   main(sys.argv)
