#!/usr/bin/python 
import argparse
import sys
import subprocess
import itertools

from Bio import SeqIO
from Bio import AlignIO

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio.Align.Applications import MuscleCommandline

import cStringIO as StringIO

from utils.threadpool import ProducerConsumer
from utils.cigar import compressCigar


def msablocks(fastaFile):
    """
    read the MSA generated by usearc/vsearch.  Each MSA ends with a consensus.
    Use this consensus sequence to identify when a MSA record ends.
    """
    accumulator= []
    for seq in SeqIO.parse(fastaFile, "fasta"):
        seq.seq = Seq(str(seq.seq).replace("-",""))
        accumulator.append(seq)
        if seq.id.lower().startswith("consensus"):
            yield accumulator
            accumulator = []
    if accumulator:
        yield accumulator


def producer(info):
    """
    A thread pool worker function.  
    - It will run muscle using stdin and stdout (no intermediate files required!
    - Build a consensus from the alignment Muscle finds
    - Drop ambigious clusters if requested
    - Build a new cigar based on the MSA alignment and the consensus that was created.
    - works on a single MSA and returns 1 consensus, with multiple cigards (1 for each sequence)
    """
    msa = info[0]
    dropAmbiguous = info[1]
    tag = info[2]
    maxhrs = info[3]
    largesize = info[4]
    newcigs = {}
    cline = ""

    if largesize != 0 and  (len(msa) - 1) >= largesize:
        if maxhrs:
            cline = MuscleCommandline(quiet = True, maxhours = maxhrs, maxiters = 2)
        else:
            cline = MuscleCommandline(quiet = True, maxiters = 2)
    else:
        if maxhrs:
            cline = MuscleCommandline(quiet = True, maxhours = maxhrs)        
        else:
            cline = MuscleCommandline(quiet = True)        

    child = subprocess.Popen(str(cline),
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True,
                             shell=(sys.platform!="win32"))
    tmp = StringIO.StringIO()   
    SeqIO.write(msa[:-1], tmp, "fasta")
    muscleout, merr = child.communicate(tmp.getvalue())
    alns = AlignIO.read(StringIO.StringIO(muscleout), "fasta")
    si = AlignInfo.SummaryInfo(alns)
    cons = si.dumb_consensus(threshold = 0, ambiguous = 'N') # take whatever the most abundant base is, that isn't a . or -
    conid =  msa[0].id.replace("*","")  + "_" + tag

    if cons.find("GATC") == -1 and dropAmbiguous:
        return []
    for s in alns:
        cigar = ""
        for i in str(s.seq):
            if i == "-" or i == '+':
                cigar += 'I' 
                #cigar += 'D' # the gap forms in the read, not the reference, therefore it should be a Deletion
            else:
                cigar += 'M'
        s.id = s.id.replace("*","")
        newcigs[s.id] = [s.id, conid, compressCigar(cigar)]
    newCon = SeqRecord(cons, id = conid, description="")
    return newCon, newcigs


def consumer(args, returndata):
    newcigs = {}
    with open(args.output1, 'w') as o:
        for r in returndata:
            if not r:
                continue
            SeqIO.write([r[0]], o, "fasta")
            newcigs.update(r[1])
    modCigars(newcigs, args.input2, args.output2)


def modCigars(newcigs, infile, outfile):
    with open(outfile, "w") as o:
        keys = set(newcigs.keys())
        for l in open(infile):
            l = l.strip().split()
            l[0] = l[0].replace("*","")
            newcigs[l[0]].append(l[-2])
            newcigs[l[0]].append(l[-1])
            print >> o, "\t".join(map(str,  newcigs[l[0]]))
            keys.remove(l[0])
        for k in keys:
            newcigs[k].append("+")
            newcigs[k].append("+")
            print >> o, "\t".join(map(str,  newcigs[k]))            


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--threads', type = int, default = 2, help = "Number of processing threads. (default: 2)" )
    parser.add_argument('-i1', '--input1',  required = True, help = "MSA output from vsearch/usearch")
    parser.add_argument('-o1', '--output1', required = True, help = " MSA consensus output")
    parser.add_argument('-i2', '--input2',  required = True, help = "userout from vsearch/usearch")
    parser.add_argument('-o2', '--output2', required = True, help = "modified userout")
    parser.add_argument('-g', '--tag', required = True, help = "what to tag consensus ID with")
    parser.add_argument('-m', '--maxhours', required = False, type = float, default = 0.0, 
                        help = "The maximum time to let an instance of muscle run (default: unlimited)")
    parser.add_argument('-n', '--largecluster', required = False, type = int, default = 0, 
                        help = "The size of a cluster that we consider too large, and will cause musle to use less sensitive parameters (default: infinite)")
    parser.add_argument("-d", "--dropAmbiguous", action = "store_true", required = False, help = "Drop clusters in which the consensus does not  contain GATC")

    args = parser.parse_args()

    worker = ProducerConsumer(args, args.threads, producer, consumer)
    worker.run( args, 
                itertools.izip( msablocks(args.input1), itertools.repeat(args.dropAmbiguous), 
                                itertools.repeat(str(args.tag)), itertools.repeat(args.maxhours) , itertools.repeat(args.largecluster)  ,))
