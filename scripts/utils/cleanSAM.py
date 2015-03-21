from Bio import SeqIO
import pysam 
import sqlite3
import cStringIO as StringIO
import os

from utils import removeFiles
from cigar import expandCigar, compressCigar
from samToBam import samToBam
from threadpool import ProducerConsumer


def producer(info):
    fileidx = str(info[0])   
    consensus = ">Consensus\n%s\n"%(str(info[2]))    
    seqs = [">%s\n%s"%(str(i), str(s)) for i, s in zip(info[3].split("\t"), info[4].split("\t"))]
    
    fastafile = StringIO.StringIO(consensus + "\n".join(seqs) )
    
    bamfile = "%s.bam"%(fileidx)
    bamidxfile = "%s.bam.bai"%(fileidx)
    with open(bamfile, "wb") as o:
        o.write(info[5])
    with open(bamidxfile, "wb") as o:
        o.write(info[6])

    refs = [ (r.id, len(r.seq) ) for r in  SeqIO.parse(fastafile, "fasta")]

    sfile = pysam.Samfile(bamfile)

    hdr = sfile.header.copy()
    hdr['SQ'] = [ {'LN': refs[0][1], 'SN': refs[0][0] }]
    
    total = 0
    indices = eval("[" + info[1] + "]")
    maxidx = max(indices) + 1
    shift = [0]*(maxidx)
    for i in xrange(len(shift)):
        shift[i] = total
        if i not in indices:
            total += 1
    outfile = pysam.Samfile("%s.trim.sam"%(fileidx), "wh", header = hdr)
    for read in sfile.fetch():
        newseq = ""
        newcig = ""
        newqual = ""
        cigar = expandCigar(read.cigarstring)
        back = 0
        pos = read.pos
        for idx, c in enumerate(cigar):
            # stop if we exceed the reference length
            if pos >= maxidx:
                break
            pos += 1
            # only process bases that we havent removed in cleanup
            # also, count any D's that we missed to shift us back in the query and qual
            if (idx + read.pos) not in indices:
                # an idx not in indices indicates that trimal has removed this particular position
                if cigar[idx].upper() == 'D':
                    back += 1
                continue

            # if the idx exists in the logs of trimAL, the base should be kept.
            newcig += cigar[idx]
            #D means we dont have a base to add.
            if cigar[idx].upper() != 'D': 
                newseq += read.query[idx - back]
                if read.qqual:
                    newqual += read.qqual[idx - back]
            else: 
                back += 1

        #if a sequence has nothing to add, why add it ?!?
        if newseq:
            read.pos = max( (read.pos - shift[read.pos]), 0)
            read.seq = newseq
            read.cigarstring = compressCigar(newcig)
            if not newqual:
                read.qual = "I"* len(newseq)
            else:
                read.qual = newqual
            outfile.write(read)
    outfile.close()
    
    samout = "%s.trim.sam"%(fileidx)
    bamout = "%s.trim.bam"%(fileidx)
    samdat = open(samout).read()
    removeFiles([bamfile, bamidxfile, samout])

    bamdat, bamidxdat = samToBam(samdat, "%s.trim"%(fileidx), buffers = True)
    return fileidx, samdat, bamdat, bamidxdat


def consumer(con, returndata):
    curs = con.cursor() 
    for dat in returndata:
        if not dat:
            continue
        sam = dat[1]
        bam = dat[2]
        bamidx = dat[3]
        ident = dat[0]
        curs.execute("""INSERT INTO trimmed_inferSAM(fileID, sam, bam, bamidx) VALUES(?,?,?,?);""", (ident, sam, sqlite3.Binary(bam), sqlite3.Binary(bamidx),))
    con.commit()


def cleanSAM(args):
    con = sqlite3.connect(args.database, check_same_thread=False)
    con.execute("""PRAGMA foreign_keys = ON;""")
    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_inferSAM(id INTEGER PRIMARY KEY, fileID INTEGER, sam TEXT, bam BLOB, bamidx BLOB, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_infersam_fileid_idx ON trimmed_inferSAM(fileID ASC);""")
    con.commit()

    curs = con.cursor()  
    curs.execute("""SELECT A.fileID, A.positions, B.sequence, C.SEQS, C.IDs, D.bam, D.bamidx 
                    FROM trimmed_logs AS A 
                         JOIN trimmed_consensus AS B ON (A.fileID = B.fileID) 
                         JOIN (SELECT fileID, group_concat(seqID, '\t') AS IDs, group_concat(sequence, '\t') AS SEQS 
                               FROM trimmed_csr GROUP BY fileID) AS C ON (C.fileID = A.fileID) 
                         JOIN inferSAM AS D ON (D.fileID = A.fileID)""")  

    rows = ( (r[0], r[1], r[2], r[3], r[4], bytearray(r[5]), bytearray(r[6]),) for r in curs)   
    worker = ProducerConsumer(args, args.threads, producer, consumer)
    worker.run(con, rows)
    con.commit()
    con.close()
