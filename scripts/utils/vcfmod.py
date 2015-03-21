#!/usr/bin/python

import vcf
import sys
import os
import argparse
import json
from collections import Counter, defaultdict
import pysam
import sqlite3
import cStringIO as StringIO

from utils.threadpool import ProducerConsumer
from utils.utils import removeFiles


def consumer(con, returndata):
    cur = con.cursor()
    for dat in returndata:
        if not dat:
            continue
        fileidx = dat[0]
        modvcf = dat[1]
        jsonstr = dat[2]
        cur.execute("""INSERT INTO trimmed_modvcf(fileID, vcf, json) VALUES (?,?,?)""", (fileidx, modvcf, jsonstr,) )
    con.commit()


def producer(info):
    fileid = str(info[0])
    vcfInput = vcf.Reader(StringIO.StringIO(info[1]))
    line = None
    try:
        line = vcfInput.next()
    except:
        return None
    if not line:
        return None

    bamfile = "%s.bam"%(fileid)
    bamidxfile = "%s.bam.bai"%(fileid)
    with open(bamfile, "wb") as o:
        o.write(info[2])
    with open(bamidxfile, "wb") as o:
        o.write(info[3])


    vcfInput = vcf.Reader(StringIO.StringIO(info[1]))
    vcfohndl = StringIO.StringIO()
    vcfOutput = vcf.Writer(vcfohndl, vcfInput)
    jsonhndl =  StringIO.StringIO()
    data = computeData(vcfInput, vcfOutput, bamfile, 0)
    json.dump(data, jsonhndl, separators=(',', ':'))
    jsonhndl.flush()
    jsonstr = jsonhndl.getvalue()
    jsonhndl.close()
    vcfohndl.flush()
    modvcf = vcfohndl.getvalue()
    vcfohndl.close()

    removeFiles([bamfile, bamidxfile])

    return info[0], modvcf, jsonstr
    

# A method to use the BAM file + the vcf to compute our data, if the vcf does not contain it, 
# or we have no idea about the specific field we need to use in the vcf file.
def computeData(inputvcf, outputvcf, samf, shift):
    samfile = pysam.Samfile(samf)
    vcff = dict([ (v.start - shift, v) for v in inputvcf])
    data = []
    for p in samfile.pileup():
        if p.pos not in vcff:
            continue
        if outputvcf:
            vcff[p.pos].POS -= shift
        if not vcff[p.pos].is_snp:
            if outputvcf:
                outputvcf.write_record(vcff[p.pos])
            continue
        ref = vcff[p.pos].REF
        alts = [str(v) for v in vcff[p.pos].ALT]
        counts = defaultdict(Counter)
        found_alts= set()
        for pread in p.pileups:
            if pread.is_del:
                continue
            base = pread.alignment.seq[pread.qpos]
            # Do we discard gaps and N?
            if base.upper() != 'N' or base != '-':
                sample = pread.alignment.qname.split("_")[-1]
                found_alts.add(base.upper())
                counts[sample].update(base.upper())
        newALTS = list(found_alts.difference([ref]))
        data.append( dict(pos=vcff[p.pos].POS, counts= dict(counts), ref = ref, alts = newALTS) )

        if outputvcf:
            vcff[p.pos].ALT = [vcf.model._Substitution(j) for j in newALTS]
            outputvcf.write_record(vcff[p.pos])
    return data


def updateVCFfiles(args):
    con = sqlite3.connect(args.database, check_same_thread=False)
    con.execute("""PRAGMA foreign_keys = ON;""")
    con.execute("""CREATE TABLE IF NOT EXISTS trimmed_modvcf(id INTEGER PRIMARY KEY, fileID INTEGER, vcf TEXT, json TEXT, FOREIGN KEY(fileID) REFERENCES files(id));""")
    con.execute("""CREATE INDEX IF NOT EXISTS trimmed_modvcf_fileid_idx ON trimmed_modvcf(fileID ASC);""")
    con.commit()
    
    cur = con.cursor()
    cur.execute("""SELECT A.fileID, A.vcf, B.bam, B.bamidx FROM trimmed_vcf as A JOIN trimmed_inferSAM AS B ON (A.fileID = B.fileID); """)
    rows = ( (r[0], r[1], bytearray(r[2]), bytearray(r[3]), ) for r in cur )

    worker = ProducerConsumer(args, args.threads, producer, consumer)
    worker.run(con, rows)
    con.commit()
    con.close()
