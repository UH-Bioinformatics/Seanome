#!/usr/bin/env python2.7
import math
import argparse
from collections import defaultdict
from scipy.special import binom


import vcf
import sqlite3
import cStringIO as StringIO

import numpy as np
import os    
import sys
import tempfile
import shutil
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt



def makeHistogram(vals, workdir, prefix = ""):
    r = (-1.0, 1.0,)
    plt.hist(vals, bins = 40, range = r)
    if not prefix:
        plt.savefig(os.path.join(workdir, "fst_histogram.png"), bbox_inches='tight')
    else:
        plt.savefig(os.path.join(workdir, "%s_fst_histogram.png"%(prefix)), bbox_inches='tight')


def computeFST(popLst):
    """
    [ pop1: [ {sample1 Q:#}, {sampleN Q:#} ] , pop2: [ {sample1 Q:#}, {sampleN Q:#} ] ]
    """
    Allele = []
    numReads = []
    for pop in popLst:
        vals = defaultdict(int)
        for samp in pop:
            for k, v in samp.iteritems():
                vals[k] += int(v)
        Allele.append(vals.values())
        numReads.append(sum(vals.values()))
    if len(numReads) != 2:
        return float("nan")
    totReads = float(sum(numReads))
    numReads1 = float(numReads[0])
    numReads2 = float(numReads[1])

    dw = [ float(Allele[0][0] * Allele[0][1] * 2),
           float(Allele[1][0] * Allele[1][1] * 2)]

    dAP = float((Allele[0][0] * Allele[1][1]) + (Allele[0][1] * Allele[1][0]))

    dblread1 = (2.0 * numReads1)
    dblread2 = (2.0 * numReads2)
    dbltotal = (2.0 * totReads)

    if dblread1 == 0 or dblread2 == 0 :
        return float("nan")

    SSW = (dw[0] / dblread1) + (dw[1] / dblread2)
    SSA = (((dAP + dw[0]) / dbltotal) - (dw[0] / dblread1)) + (((dAP + dw[1]) / dbltotal) - (dw[1] / dblread2))
    dfW = (numReads1 - 1.0) + (numReads2 - 1.0)
    dfA = 1.0
    MSW = SSW / dfW
    MSA = SSA / dfA
    Nc = totReads - (( (numReads1**2) + (numReads2**2) ) / totReads)
    s2A = (MSA - MSW) / Nc
    if dAP != 0:
        return s2A / (s2A + MSW)
    return 0


def processVCF(fname, stream, con):
    vcffile = vcf.Reader(fsock = stream)
    retdat = []
    samplemap = dict([ (i, e,) for e, i in enumerate(vcffile.samples)])
    samples = list(vcffile.samples)
    for i in  con.execute("""SELECT sampleid, species, A.id FROM groups AS A JOIN files AS B ON (B.id = A.fileID) WHERE B.name = ?;""", (fname,)):
        samples[ samplemap[i[0]] ] = (samples[ samplemap[i[0]] ], i[1], i[2])
    
    
    for row in vcffile:
        dat = [fname, str(row.CHROM), str(row.POS), str(row.REF), ":".join([str(q) for q in row.ALT]), ":".join(["1"]*len(row.samples))]
        process = True

        for sample in row.samples:
            data = sample.data
            if data.RD == None or data.AD == None:
                ref = 0
                alt = 0
                #process = False
                #break
            else:
                ref = str(data.RD)
                alt = str(data.AD)
            #total = data.SDP
            dat.append("%s:%s"%(ref, alt))
        if process:
            retdat.append( ("\t".join(dat), samples,))
    return retdat, samples


def main():
    parser = argparse.ArgumentParser(description='Extracts the vcf files after vcf_generator.py')
    parser.add_argument('-d', '--database', help = 'Specified Seanome Database', required = True)
    parser.add_argument('-w', '--work', help = 'Work direcory', required = True)
    parser.add_argument('-p', '--prefix', help = 'output prefix', required = False, default = "")
    args = parser.parse_args()
    
    con = sqlite3.connect(args.database, check_same_thread = False)
    cur = con.cursor()
    cur.execute(""" SELECT B.name, A.vcf, B.id FROM trimmed_vcf AS A JOIN files AS B ON (A.fileID = B.id); """)

    alldat = []
    sampleMap = []
    fileids = []
    for dat in cur:
        if not dat:
            continue
        ret, mapping = processVCF(dat[0], StringIO.StringIO(dat[1]), con)
        alldat.extend( [ (r[0], r[1], dat[2],) for r in ret])
        #alldat.extend(ret)
        #sampleMap.append(mapping)
        #fileids.append(dat[2])

    fsts = []
    fstsprint = []
    for l, idLst, fid in alldat:
    #for l, idLst, fid in zip(alldat, sampleMap, fileids):
        l = l.strip().split("\t")
        info = l[:6]
        #info[5] # contains a colon delimited string that determines how many in the vals section belong to a given population
        vals = [ map(int, q.split(":")) for q in l[6:] ]

        ref = info[3]
        alts = info[4]
        pops = []
        idx = 0
        for d in info[5].split(":"):
            d = int(d)
            for i in xrange(d):
                parts = vals[idx]
                idx += 1
                pops.append([{ref : int(parts[0]), alts : int(parts[1])}] )

        for x in range(len(pops)):
            for y in range(x+1, len(pops)):
                result = computeFST([pops[x], pops[y]])
                if not math.isnan(result):
                    fsts.append(result)
                    #print str(fid)
                    #print idLst[x]
                    #print idLst[y]
                    con.execute("""INSERT INTO fst(fileID, groupA, groupB, pos, value) VALUES(?,?,?,?,?)""", (str(fid), idLst[x][2], idLst[y][2], info[2], result,) )
    con.commit()
    con.close()
    if fsts:
        makeHistogram(fsts, args.work, prefix = args.prefix)    
    else:
        print >> sys.stderr , "No data to generate histogram from"


if __name__ == "__main__":
    #try:
    main()
    #except:
    #    print >> sys.stderr, "An error occurred while trying to create the fst histogram"
    shutil.rmtree(os.environ['MPLCONFIGDIR'])
#        print "%s\t%s\t%s\t%s\t%s"%(info[0], info[1], info[2], str(computeFst(pops)), str(computeFST(pops)) )
