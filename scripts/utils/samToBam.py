import subprocess
import os
from utils import removeFiles


def samToBam(samdat, fprefix, buffers = True):
    if not buffer:
        return samToBamNoBuffer(samdat, fprefix)
    else:
        return samToBamBuffers(samdat, fprefix)


def samToBamNoBuffer(samdat, fprefix):
    bamout = "%s.bam"%(fprefix)
    cline = """samtools view -bS -"""
    child = subprocess.Popen(str(cline),
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=(sys.platform!="win32"),
                             close_fds = True)
    sout, serr = child.communicate(samdat)

    cline = """samtools sort - %s"""%(fprefix)
    child2 = subprocess.Popen(str(cline),
                              stdin = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              stdout = subprocess.PIPE,
                              shell = (sys.platform!="win32"),
                              close_fds = True)
    bamdat, berr = child2.communicate(sout)
    os.system("""samtools index %s > /dev/null 2> /dev/null"""%(bamout))
    return bamout, "%s.bai"%(bamout)


def samToBamBuffers(samdat, fprefix):
    cline = """samtools view -bS -"""
    child = subprocess.Popen(str(cline),
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=(sys.platform!="win32"),
                             close_fds = True)
    sout, serr = child.communicate(samdat)
    bamout = "%s.tmp.bam"%(fprefix)
    cline = """samtools sort - -o %s.tmp"""%(fprefix)
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
    removeFiles([bamout, "%s.tmp.bam.bai"%(fileidx)])

    return bamdat, bamidxdat
