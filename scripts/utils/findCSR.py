import re
import subprocess
import sys
import sqlite3
from Bio import AlignIO
import tempfile
from Bio import SeqIO
import shutil

HMM_CMD="""nhmmer --cpu 1 --qformat fasta - % """
hmmer_hit_re = re.compile(r'score\s*bias\s*Evalue\s*hmmfrom\s*hmm\s*to\s*alifrom')
hmmer_query_length = re.compile("Query:\s*ali\s*\[M=(\d+)")


def find_csr(args):
   con = sqlite3.connect(args.database, check_same_thread=False)
   cur = con.cursor()
   cur.execute("""SELECT fileID, group_concat(seqID, '\t') AS IDs, group_concat(sequence, '\t') AS SEQS 
                  FROM  csr 
                  GROUP BY fileID;""")
   genomeseq = SeqIO.read(args.genome, "fasta")
   tmpDir = tempfile.mkdtemp(prefix="seanome_", dir=".")
   nhmmer_dir = os.path.join(tmpDir, "nhmmer_dir")
   os.mkdir(nhmmer_dir)
   worker = ProducerConsumer(args, args.threads, find_csrProducer, find_csrConsumer)
   worker.run(con, ( (c, args.genome, genomeseq, args.min_csr_len, args.min_csr_sim) for c in cur) )
   con.commit()
   con.close()
   # DLS cleanup after mkdtemp.  It is the Users responsibility!  
   try:
      shutil.rmtree(tmpDir, onerror = remove_readonly)
   except:
      pass


def getSubSeqFromFasta(seq, start, end, reverse=False):
    ident = "%s_%s_%s%s" % (seq.id, start + 1, end, "_R" if(reverse) else "" )
    if reverse:
        seqo = seq[start:end].reverse_complement().seq
    else:
        seqo = seq[start:end].seq
    return ident, seqo


def find_csrProducer(payload):
    info, genome, genomeseq, min_csr_len, min_csr_sim = payload
    seqs = "\n".join([">%s\n%s"%(str(i), str(s)) for i, s in zip(info[1].split("\t"), info[2].split("\t"))])
    cline = HMM_CMD%(genome)
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


def find_csrConsumer(con, returndat):
   cur2 = con.cursor()
   for dat in returndat:
       if not dat:
           continue
       print dat
       for q in dat[1]:
           print str(q)
           if q.id == dat[2]:
               cur2.execute("""INSERT INTO csr(fileID, seqID, sequence) VALUES(?, ?, ?)""", (dat[0], q.id, str(q.seq), ) )
           else:
               cur2.execute("""UPDATE csr SET sequence = ?  WHERE seqID = ?;""", (str(q.seq), q.id,) )
   con.commit()
