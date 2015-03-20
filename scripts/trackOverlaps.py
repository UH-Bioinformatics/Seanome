#!/usr/bin/python
import re
import sys
from collections import defaultdict
from Bio import SeqIO
import argparse

from utils.cigar import expandCigar, compressCigar


def modifyCig(cig):
    """
    convert all leading and trailing I and D's into M.  
    the parseClusteringOut has does not distinguish and 
    will blindly fill all + - in with some type of base.
    For now, we will just make all Is and Ds int Ms 
    """
    cig = list(cig)
    for idx in xrange(len(cig)):
        if cig[idx] == 'M':
            break
        cig[idx] = 'M'
    #print "".join(cig)
    cig = cig[::-1]
    for idx in xrange(len(cig)):
        if cig[idx] == 'M':
            break
        cig[idx] = 'M'
    cig = cig[::-1]
    #print "".join(cig)
    return "".join(cig)    


def updateCigar(child, childstrand, newcentroid, centroidstrand, iterNum, newCigar, seqConsCigar, outputfile, singlestr = ""):
    process = [ [iterNum, child, s, idx, False, newCigar, False if(childstrand !=   '-') else True, newcentroid,] for idx, s in enumerate(seqConsCigar[iterNum][child])]
    while process:
        iNum, ident, info, idx, isNewCentroid, newCigar, isreverse, newparent = process.pop(0)
        if info[0] == newparent:
            continue
        myCigar = info[1]
        orignewcig = newCigar 
        if isreverse:
            myCigar = info[1][::-1]
        posCigar = 0
        posRef = 0
        posRead = 0
        newCig = ""
        while posRef < len(newCigar) and posRead < len(myCigar):
            if newCigar[posRef] == 'I':
                newCig += 'I'
                posRef += 1
            elif newCigar[posRef] == 'D':
                # if we come across a 'D' in the new cigar, we need to keep it
                newCig += 'D'
                posRef += 1
            elif myCigar[posRead] == 'I':
                newCig += 'I'
                posRead += 1
                posRef += 1
            else: 
                newCig += 'M'
                posRef += 1
                posRead += 1
        while posRef < len(newCigar):
            newCig += "I"
            posRef += 1
        if isreverse and not isNewCentroid:
            if seqConsCigar[iNum][ident][idx][2] == '+':
                seqConsCigar[iNum][ident][idx][2] = '-'
                info[2] = '-'
            else:
                seqConsCigar[iNum][ident][idx][2] = '+'
                info[2] = '+'

            seqConsCigar[iNum][ident][idx][3] = '+'
            info[3] = '+'
        seqConsCigar[iNum][ident][idx][1] = newCig

        print >> outputfile, "%s\t%s\t%s\t%s\t%s\t*"%(info[0], newparent, compressCigar(newCig), info[2], info[3])
        info[1] = newCig


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--input1',  required = True, help = "output of muso _1" )
    parser.add_argument('-i2', '--input2',  required = True, help = "output of muso _2" )
    parser.add_argument('-o', '--output',  required = True, help = "output file")
    args = parser.parse_args()
    seqConsCigar = {}

    iterNum = 0
    fileName = args.input1
    seqConsCigar[iterNum] = defaultdict(list)
    for line in open(fileName, 'r'):
        # each line is:: s1     s2     cigar       strand_s1       strand_s2
        data = line.strip().split()
        seqConsCigar[iterNum][data[1]].append([data[0], expandCigar(data[2]), data[3], data[4]])

    iterNum = 1
    fileName = args.input2

    with open(args.output, "w") as o:
        seqConsCigar[iterNum] = defaultdict(list)
        for line in open(fileName, 'r'):
            # each line is:: s1     s2     cigar       strand_s1       strand_s2
            data = line.strip().split()
            if seqConsCigar[iterNum-1].has_key(data[0]): # we are intersted in the sequence that hit against the new centroid
                updateCigar(data[0], data[3], data[1], data[4],   iterNum-1, expandCigar(data[2]), seqConsCigar, o)         
