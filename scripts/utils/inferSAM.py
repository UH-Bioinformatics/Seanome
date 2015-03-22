import os
import copy

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
import re
import sys
import cStringIO as StringIO
import subprocess
from collections import Counter

from sqlitedb import QUERY_CSR_AS_SEQS
from utils import removeFiles, CONSENSUS_NAME
from threadpool import ProducerConsumer


def inferSAMProducer(info):
    name = str(info[0])
    fastseqs = info[1] + info[2]
    samdir = info[3]
    sb = SAM_BUILDER(name, fastseqs, samdir)
    return sb.run()


def inferSAMConsumer(con, returndata):
    curs = con.cursor()
    for dat in returndata:
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


def inferSAM(args, con):
    con = sqlite3.connect(args.database, check_same_thread=False)
    curs = con.cursor()
    curs.execute("""SELECT A.fileID, B.sequence, A.IDs, A.SEQS 
                    FROM ( %(csr_as_seqeuences)s ) AS A 
                    JOIN consensus AS B ON (A.fileID = B.fileID);"""%sict(csr_as_seqeuences = QUERY_CSR_AS_SEQS) )
    rows = (  (r[0], ">%(seqID)s\n%(seq)s\n"%dict(seqID = CONSENSUS_NAME, seq[1]), "\n".join([">%s\n%s"%(i,s,)        for i, s in itertools.izip(r[2].split("\t") , r[3].split("\t") ) ] ) , sd,) for r, sd in 
              itertools.izip(curs, itertools.repeat(str(args.split_sam_dir)) )    )

    worker = ProducerConsumer(args, args.threads, inferSAMProducer, inferSAMConsumer)
    worker.run(con, rows)



class SAM_BUILDER():
    cigar_re = re.compile(r'([0-9]+)([M=XID])')
    soft_clip = re.compile(r'(^[0-9]+)([S])')


    def __init__(self, name, fastseqs, samdir):
        self.name = name
        self.data = fastseqs
        self.samdir = samdir

    def getShiftedPosRef(self, seq, startInRef):
        pos = 0
        nbMatches = 0
        if startInRef == 0:
            return pos
        #print startInRef
        #print len(seq)
        while nbMatches < startInRef:
            if seq[pos] != "-":
                nbMatches += 1
                pos += 1
            else:
                pos += 1
        return pos


    def getShiftedPosRead(self, cigar, start):
        pos = 0
        shiftedReadStart=start
        nbMatches = 0
        while nbMatches < start:
            if cigar[pos] == "I" :
                shiftedReadStart += 1
            else:
                nbMatches+=1
            pos += 1
        return shiftedReadStart


    def expandCigar(self, cigar):
        cigstr = ""
        for m in self.cigar_re.finditer(cigar):
            cigstr +=  m.group(2) * int(m.group(1))
        return cigstr


    def compressCigar(self,cigar):
        cigar = cigar.rstrip("D")
        c = None
        cnt = 0
        ciglst = []
        for b in cigar:
            if b != c:
                if c!= None:
                    ciglst.append("%s%s"%(cnt, c)) 
                cnt = 1
                c = b
            else:
                cnt += 1
        if c!= None:
            ciglst.append("%s%s"%(cnt, c))     
        return "".join(ciglst) 


    def __getTag__(self, sinfo, tag):
        try:
            return sinfo.opt(tag)
        except KeyError:
            return None


    def generateSamInfo(self,  refnamestr, sinfo, seq, cig, startpos, refname, qual, readgroupID):
        a = pysam.AlignedRead()
        a.rname = refname
        a.qname = "%s_%s"%(sinfo.qname, refnamestr)
        a.seq= seq
        a.flag = 0x00 #sinfo.flag & 16
        a.pos =  startpos
        a.mapq = sinfo.mapq
        a.cigarstring =  cig
        a.rnext = -1
        a.pnext= -1 #a.pos
        #a.isize = sinfo.isize
        a.tlen = 0
        a.qual = qual.strip().replace(" ","")
        tags=  []
        # DLS: force to use our regroup information instead of what was previously filled in..
        tags.append( ("RG", readgroupID,) )
            
        tmp = self.__getTag__(sinfo, "X0")
        if tmp:
            tags.append( ("X0", tmp,) )
        tmp = self.__getTag__(sinfo, "AS")
        if tmp:
            tags.append( ("AS", tmp,) )
        tmp = self.__getTag__(sinfo, "XS")
        if tmp:
            tags.append( ("XS", tmp,) )
        tmp = self.__getTag__(sinfo, "YS")
        if tmp:
            tags.append( ("YS", tmp,) )
        tmp = self.__getTag__(sinfo, "YT")
        if tmp:
            tags.append( ("YT", tmp,) )
        a.tags = tags
        #a.tags = sinfo.tags
        return a


    def parseMSA(self, ffile):
        msaRefs = dict()
        for seq in  SeqIO.parse(ffile, 'fasta'):
            if not seq.id.startswith("Con"):
                nameInfo = seq.id.split("_")
                msaRefs[nameInfo[0]]=[int(nameInfo[1])-1, int(nameInfo[2])-1, str(seq.seq)]
                if nameInfo[-1] == 'R': # changed to look at the last field for a capital R instead of the header length..
                #if len(nameInfo)==4:
                    msaRefs[nameInfo[0]].append("Reversed")
        return msaRefs

    
    def computeStartPositions(self, read, refStart, refEnd, isReverse):
        if not isReverse: # this means ref is forward
            #if read starts Before msaRef
            if read.pos < refStart:
                startInRef = 0
                startInRead = refStart - read.pos
            else:
                startInRead = 0
                startInRef = read.pos - refStart
        else:
            endbase = read.aend - 1
            if endbase < refEnd:
                startInRef  =  refEnd - endbase
                startInRead =  0
            else:
                startInRef  =  0 
                startInRead =  endbase - refEnd
        return startInRef, startInRead


    def processSingleRef(self, samfile, name, refStart, refEnd, refSeq, isReversed, combinedOut, header, readgroupID):
        myReads = samfile.fetch(name, refStart, refEnd)

        for read in (q for q in myReads if not q.is_unmapped):
            startInRef, startInRead = self.computeStartPositions(read, refStart, refEnd, isReversed)
            #print refSeq
            #print read.qname
            shiftedStartRef = self.getShiftedPosRef(refSeq, startInRef)
            expandedCigar  = self.expandCigar(read.cigarstring)
            shiftedStartRead = self.getShiftedPosRead(expandedCigar, startInRead)
            
            quality = ""
            # DLS 20150303  Addition of counting D to subtract from where we cut a sequence.
            if not isReversed:
                cutrgn = expandedCigar[:shiftedStartRead]
                readCigar = expandedCigar[shiftedStartRead:]
                backped = cutrgn.count("D")
                readString = read.query[(shiftedStartRead - backped):]
                if read.qqual:
                    quality = read.qqual[(shiftedStartRead - backped):]
            else:
                revcigar = expandedCigar[::-1]
                cutrgn = revcigar[:shiftedStartRead]
                readCigar = revcigar[shiftedStartRead:]
                backped = cutrgn.count("D")
                readString = str(Seq(read.query).reverse_complement()[(shiftedStartRead - backped):])
                if read.qqual:
                    quality = read.qqual[::-1][(shiftedStartRead - backped):]

            refString = refSeq[shiftedStartRef:]          
            #DEBUGGING
            # print "read is %s" % read.qname
            # print "orignal Cigar is %s " % read.cigarstring
            # print "original position %s " % read.pos
            # print "original read str %s " % read.query
            
            # print "startInRef is %s" % (startInRef)
            # print "shiftedStartRef is %s"  % shiftedStartRef

            # print "startInRead is %s"  % startInRead
            # print "shiftedStartRead is %s" % (shiftedStartRead)

            # print readCigar
            # print refString
            # print readString
            # inspecting columns and inserting gaps or deleting Inserts from read as necessary
            posRef = 0
            posRead = 0
            posCigar = 0
            newReadString = ""
            newCigarString = ""
            tempString = "" # Just visualization purposes for now
            newQualString = ""

            while posRef < len(refString) and posRead < len(readString):
                #print refString[:posRef+1]
                #print readString[:posRead+1]
                #print readCigar[:posCigar+1]
                if refString[posRef] == "-":
                    if readCigar[posCigar] == "I":
                        newReadString += readString[posRead]
                        tempString += readString[posRead]
                        if read.qqual:
                            newQualString += quality[posRead]
                        newCigarString += "M"
                        posRef += 1
                        posRead += 1
                        posCigar += 1
                    elif readCigar[posCigar] == "D":
                        newCigarString += readCigar[posCigar]
                        posRef += 1
                        #posCigar += 1
                    else:
                        tempString += "-"
                        newCigarString += "D"
                        posRef += 1                    
                elif readCigar[posCigar] == "I": # we know ref has a non-gap
                    posRead += 1
                    posCigar += 1
                    #print "I encountered"
                elif readCigar[posCigar] == "D":
                    tempString += "-"
                    newCigarString += readCigar[posCigar]
                    posRef += 1
                    posCigar += 1
                else:
                    newReadString += readString[posRead]
                    if read.qqual:
                        newQualString += quality[posRead]
                    tempString += readString[posRead]
                    newCigarString += readCigar[posCigar]
                    posRef += 1
                    posRead += 1
                    posCigar += 1
                #print tempString
                #print ""
                #raw_input()

            compactCigar = self.compressCigar(newCigarString)
            samData = self.generateSamInfo(name, read, newReadString, compactCigar, shiftedStartRef, 0, newQualString, readgroupID)

            #DEBUGGING
            #print "refStr : %s" % refString
            #print "tempStr: %s" % tempString
            #print "cigrStr: %s" % newCigarString
            #print "cmpCigr: %s\n" % compactCigar
            
            combinedOut.write(samData) 


    def generateReadGroups(self, groups):
        RG_template = { 'ID': '',
                         'LB': '', 
                         'SM': '',
                         'PL': 'ILLUMINA'} # Platform/technology used to produce the reads. 
        readgroups = []
        idtomap = {}
        for idx, g in enumerate(groups):
            rgrp = RG_template.copy()
            rgrp['ID'] = "GROUP-"+str(idx + 1)
            rgrp['LB'] = g
            rgrp['SM'] = g
            readgroups.append(rgrp)
            idtomap[g] = "GROUP-"+str(idx + 1)
        return readgroups, idtomap


    def updateConsensus(self, bamfile):
        samfile = pysam.Samfile(bamfile)
        colBases = []
        for pileupcolumn in samfile.pileup(CONSENSUS_NAME):
            bases =[]
            for pup in pileupcolumn.pileups:
                bases.append(pup.alignment.seq[pup.qpos])
            colBases.append(bases)
        newSeq=""
        for pos in colBases:
            values = sorted(Counter(pos).items(), key=lambda x: x[1], reverse=True)
            added = False
            for val in values:
                if val[0] == "N" or val[0] == "-":
                    continue
                else:
                    newSeq += val[0]
                    added = True
                    break
            #DLS we need to add a base.. IF somehow we end up with a column with no 
            #    majority, other than N and -, we need to put N...
            if not added:
                newSeq += "N"
        return newSeq


    def run(self):
        samout = "%s.sam"%(self.name)
        bamout = "%s.bam"%(self.name)

        msaInfo = self.parseMSA(StringIO.StringIO(self.data))

        refs = [ (r.id, len(r.seq) ) for r in  SeqIO.parse(StringIO.StringIO(self.data), "fasta")]

        readgroups, refIDToReadgroup = self.generateReadGroups(msaInfo.keys())
        header = dict(HD = dict(VN = '1.0'), SQ = [ {'LN': refs[0][1], 'SN': refs[0][0] }], RG = readgroups)
        outfile = pysam.Samfile(samout, "wh", header = header)       
        for refName, refVals in msaInfo.iteritems():
            isReversed = True
            if len(refVals) == 3:
                isReversed = False
            samid = os.path.join(self.samdir,"%s.bam" % refName )
            samfile = pysam.Samfile(samid, "rb" )        
            self.processSingleRef(samfile, refName, refVals[0], refVals[1], refVals[2], isReversed, 
                                  outfile, header, refIDToReadgroup[refName])
            samfile.close()
        outfile.close()

        # need it in memory anyways...
        samdat = open(samout, "rb").read()      
        cline = """samtools view -bS -"""
        child = subprocess.Popen(str(cline),
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=(sys.platform!="win32"),
                                 close_fds = True)
        sout, serr = child.communicate(samdat)

        cline = """samtools sort - -o %s"""%(self.name)
        child2 = subprocess.Popen(str(cline),
                                  stdin = subprocess.PIPE,
                                  stderr = subprocess.PIPE,
                                  stdout = subprocess.PIPE,
                                  shell = (sys.platform!="win32"),
                                  close_fds = True)
        bamdat, berr = child2.communicate(sout)
        with open(bamout, "wb") as o:
            o.write(bamdat)
        
        os.system("""samtools index %s > /dev/null 2> /dev/null"""%(bamout))
        bamidxdat = open("%s.bai"%(bamout), "rb").read()

        newcon = self.updateConsensus(bamout)
        removeFiles([samout, bamout, "%s.bai"%(bamout)])
        return samdat, bamdat, bamidxdat, self.name, refIDToReadgroup.values(), newcon
