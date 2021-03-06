import pysam
from collections import Counter


def updateConsensus(bamfile):
    samfile = pysam.Samfile(bamfile)
    colBases = []
    for pileupcolumn in samfile.pileup():
        bases = []
        for pup in pileupcolumn.pileups:
            try:
                bases.append(pup.alignment.seq[pup.qpos])
            except:
                if pup.query_position:
                    bases.append(pup.alignment.seq[pup.query_position])
                #else:
                #    bases.append("-")
        colBases.append(bases)
    newSeq=""
    for pos in colBases:
        values = sorted(Counter(pos).items(), key = lambda x: x[1], reverse = True)
        added = False
        for val in values:
            if val[0] == "N" or val[0] == "-":
                continue
            else:
                newSeq += val[0]
                added = True
                break
        if not added:
            newSeq += "N"
    return newSeq
