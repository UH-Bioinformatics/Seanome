import sqlite3
from tempfile import mkstemp
import cStringIO as StringIO
from Bio import SeqIO
from utils import removeFiles
import os
import subprocess 
import sys

USRCH= """usearch -usearch_local %(input)s -db %(database)s -id 0.80 -strand both -userfields query+target+id+alnlen+bits+qstrand+tstrand+qlo+qhi+tlo+thi+qrow+trow -userout %(output)s -threads %(threads)i -evalue 1e-5 %(tail)s"""



def updateSearchOutput(args, fastaOne, fastaTwo, searchdata):
    dataQ = {}
    dataT = {}
    pos = 0
    searchdata.sort( key = lambda x: (x[0], int(x[3]),), reverse = True )
    for seq in SeqIO.parse(fastaOne, "fasta"):
        end = pos + len(seq.seq)
        dataQ[seq.id] = pos
        pos = end
    pos = 0
    for seq in SeqIO.parse(fastaTwo, "fasta"):
        end = pos + len(seq.seq)
        dataT[seq.id] =pos
        pos = end
    seen = set()

    data = StringIO.StringIO()
    for l in searchdata:
        if l[0] in seen:
            continue
        seen.add(l[0])
        l[7]  = str(int(l[7]) + dataQ[l[0]])
        l[8]  = str(int(l[8]) + dataQ[l[0]])
        l[9]  = str(int(l[9]) + dataT[l[1]])
        l[10] = str(int(l[10]) + dataT[l[1]])
        print >> data, "\t".join(l)
    return data.getvalue()


def find_shared_regions(args):
    tmpname = mkstemp(dir = ".")
    os.close(tmpname[0])
    tmpname = tmpname[1]
    qaction = ""
    if args.quiet:
        qaction = " > /dev/null 2> /dev/null "
    cline = USRCH%dict(input=args.input1, database=args.input2, output= tmpname, threads = args.threads, tail = qaction )
    print cline
    child = subprocess.Popen(str(cline), shell=(sys.platform!="win32") )
    child.wait()
    data = [l.strip().split() for l in open(tmpname)]
    removeFiles([tmpname])
    return data


def find_seed_csr(args, con):   
    inFile = updateSearchOutput( args,args.input1, args.input2, find_shared_regions(args)) 

    filemap = dict()
    outFileNum=0
    c = con.cursor()
    for line in inFile.splitlines():
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
                seqid1 = "%s_%s_%s%s"%(args.name1, data[7], data[8] , "_R" if(data[5] == "-") else "" )
                seqid2 = "%s_%s_%s%s"%(args.name2, data[9], data[10], "_R" if(data[6] == "-") else "" )
            c.execute("""INSERT INTO csr(fileID, seqID, sequence, seed) values(?,?,?,?)""", (fileID, seqid1, data[11], 1) )
            c.execute("""INSERT INTO csr(fileID, seqID, sequence, seed) values(?,?,?,?)""", (fileID, seqid2, data[12], 1) )
            outFileNum += 1
    con.commit()
