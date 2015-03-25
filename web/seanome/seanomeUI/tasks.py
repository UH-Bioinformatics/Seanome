from __future__ import absolute_import
import os
import sys
import time
import threading
import subprocess
from collections import defaultdict
from celery import shared_task

from celery.utils.log import get_task_logger
logger = get_task_logger(__name__)

from .models import *
from django.conf import settings
from Bio import SeqIO
from cStringIO import StringIO
from django.db.models import Q
from django.db import connection, close_old_connections
from redis import Redis
from redis_semaphore import Semaphore

import json
import datetime


@shared_task
def cleanupSeanome(removeOlder):
    try:
        removeOlder = int(removeOlder)
    except:
        removeOlder = 10
    thresh = datetime.datetime.utcnow() - datetime.timedelta(days=removeOlder)

    for job_object in job.objects.filter(Q(modifydate__lte = thresh)).exclude(Q(stage__lt = 4)):
        # remove job site
        pth = os.path.join(settings.JOBDIR, str(job_object.id))
        os.system("rm -r %s"%(pth))
    job.objects.filter(Q(modifydate__lte = thresh)).exclude(Q(stage__lt = 4)).delete()
    return True


@shared_task
def seanomeBad(directory):
    """
    Failed validation of IDs..
    """
    os.system("""rm -rf '%s'"""%(directory))
    return True


def reorder(row):
    reads = []
    asm = []
    for rf in row:
        for l in open(rf):
            if l.strip().startswith(">"):
                asm.append(rf)
                break
            elif l.strip().startswith("@"):
                reads.append(rf)
                break
    return reads + asm


def makeStage1Dirs(workdir, skipmaskrepeat, action):
    # os.system("mkdir %s"%(os.path.join(workdir, "00_inputs")))
    # if action != '2':
    #     os.system("mkdir %s"%(os.path.join(workdir, "01_pear")))
    #     os.system("mkdir %s"%(os.path.join(workdir, "02_sort")))
    #     if action == '1':
    #         os.system("mkdir %s"%(os.path.join(workdir, "03_cap")))
    #     os.system("mkdir %s"%(os.path.join(workdir, "03_clustered")))
    #     if skipmaskrepeat != '1':
    #         os.system("mkdir %s"%(os.path.join(workdir, "04_dusted")))
    #         os.system("mkdir %s"%(os.path.join(workdir, "05_dedup")))
    #         os.system("mkdir %s"%(os.path.join(workdir, "06_repeats")))
    #         os.system("mkdir %s"%(os.path.join(workdir, "07_usearch")))
    # else:
    #     os.system("mkdir %s"%(os.path.join(workdir, "01_combined")))
    #     os.system("mkdir %s"%(os.path.join(workdir, "02_coverage")))
    #     if skipmaskrepeat != '1':
    #         os.system("mkdir %s"%(os.path.join(workdir, "03_dusted")))
    #         os.system("mkdir %s"%(os.path.join(workdir, "04_dedup")))
    #         os.system("mkdir %s"%(os.path.join(workdir, "05_repeats")))
    #         os.system("mkdir %s"%(os.path.join(workdir, "06_usearch")))
    # os.system("mkdir %s"%(os.path.join(workdir, "08_clean")))
    # os.system("mkdir %s"%(os.path.join(workdir, "09_stats_s1")))
    os.system("""chmod -R 755 %s"""%(workdir))


 def makeStage2Dirs(workdir, skipmaskrepeat, action):
     pass
#     os.system("mkdir %s"%(os.path.join(workdir, "10_uclust_to_sam")))
#     os.system("mkdir %s"%(os.path.join(workdir, "12_maskrepeats")))
#     if skipmaskrepeat != '1':
#         os.system("mkdir %s"%(os.path.join(workdir, "11_bowtie2")))
#     os.system("""chmod -R 755 %s"""%(workdir))


class processSubtask(threading.Thread):
    def __init__(self, command, semA, semB):
        super(processSubtask, self).__init__()
        self.command = map(str, command)
        self.semA = semA
        self.semB = semB
        self.start()

    def run(self):
        task = subprocess.Popen(self.command)
        task.wait()
        self.semA.release()
        self.semB.release()


@shared_task
def seanomeGO(workdir, jid, idmap):
    jobj = job.objects.get(pk = jid)

    jobj.stage = STAGES_REVERSE['running']
    sparams = json.loads(jobj.scriptparams)
    jobj.save()

    os.system("""chmod -R 755 %s"""%(workdir))
    setupUserEnv(workdir)
    close_old_connections()

    #skipmasking = "0" if(sparams.get('repeatmask', 0) == 1) else "1" 
    
    # just in case something very odd occurs and the script fails to make the directories...
    makeStage1Dirs(workdir, skipmasking, action)
    
    tasks = []
    libs = {}
    coverageFiles =[]
    for k, row in idmap.iteritems():
        infile_F = os.path.join(workdir, row[0])
        infile_R = os.path.join(workdir, row[1])        
        try:
            infile_Ftmp = os.path.join(workdir, row[0].replace(" ","_"))
            infile_Rtmp = os.path.join(workdir, row[1].replace(" ","_"))
            os.rename(infile_F, infile_Ftmp)
            os.rename(infile_R, infile_Rtmp)
            infile_F = infile_Ftmp
            infile_R = infile_Rtmp
            libs[k] = [infile_F, infile_R]
            # We need this info for later..
            coverageFiles.append("%s"%(os.path.join(workdir, k, "%s.coverage.json"%(k)) ) )
        except:
            pass

   with open(os.path.join(workdir, "config.yaml"), "w") as ymout:
       yaml.dump(dict(libraries=libs, coverage=coverageFiles), ymout)
   # This will generate --ALL-- scripts that are requied to run the seanome pipeline.  instead of using the primary script, we will utilize the intermediate scripts
   os.system("SeanomeWrapper.py -c %s -d %s -t %s -j 1 multiple"%( os.path.join(workdir, "config.yaml") , os.path.join(workdir, "seanome.db3") , settings.WORK_THREADS))

   # used in an attempt to limit how many threads of processing a task spins up at any given time.
   rdis_semaphore = Semaphore(Redis(), count = settings.REDIS_THREADS, namespace = 'SeanomeWorkers')
   rdis_semaphore.acquire()       
   if not jobj.interrupt:
       os.system("""bash primary_script.sh >> log.out 2>> log.err""")
       jobj = job.objects.get(pk = jid)
       jobj.stage = STAGES_REVERSE['done-s2']
       jobj.save()       
   else:
       os.system("""bash child_runner_script_1.sh  >> log.out 2>> log.err """)

       jobj = job.objects.get(pk = jid)
       jobj.stage = STAGES_REVERSE['done-s1']
       jobj.save()

   rdis_semaphore.release()
   return True


@shared_task
def seanomeGOstage2(workdir, jid, cutoffs):
    # DLS -- need to also add the same option that seanome.sh has, so that we can skip
    # the repeatmasking steps.
    jobj = job.objects.get(pk = jid)
    sparams = json.loads(jobj.scriptparams)
    jobj.stage = STAGES_REVERSE['running']
    jobj.save()
    setupUserEnv(workdir)
    close_old_connections()

    skipmasking = "0" if(sparams.get('repeatmask', 0) == 1) else "1" 

    makeStage2Dirs(workdir, skipmasking, action)

    # used in an attempt to limit how many threads of processing a task spins up at any given time.
    rdis_semaphore = Semaphore(Redis(), count = settings.REDIS_THREADS, namespace = 'SeanomeWorkers')

    for l in open("child_runner_script_2.sh"):
        if l.startswith("bash "):
            rdis_semaphore.acquire()
            l = l.strip().split()
            key = l[1].split("_")[0]
            v = cutoffs[key]
            v.sort()                     
            l = " ".join(l[:2])
            if len(v) == 2:
                os.system("""%s %s %s  >> log.out 2>> log.err """%(l, v[0], v[1]))
            else:
                os.system("""%s %s %s  >> log.out 2>> log.err """%(l, 3, 200))
        rdis_semaphore.release()

    rdis_semaphore.acquire()
    os.system("""bash child_runner_script_3.sh  >> log.out 2>> log.err """)
    rdis_semaphore.release()

    jobj = job.objects.get(pk = jid)
    jobj.stage = STAGES_REVERSE['done-s2']
    jobj.save()
    return True


def setupUserEnv(workdir = None):
    if workdir:
        os.chdir(workdir)
    os.environ['PATH'] += os.pathsep + settings.USER_BIN


def createdir(path):
    try:
        os.makedirs(path)
    except OSError:
        pass


# def prepCSROutput(parent, idA, idB):
#     path = os.path.join(parent, "%s_%s"%(idA, idB))
#     createdir(path)
#     return path


# def generateOrdering(workdir, infodir, idents, params, searchmode = 1):
#     outdir = os.path.join(workdir, settings.ORDER_DIR)
#     createdir(outdir)
#     if searchmode == SEARCH_MODE_DICT['exhaustive']:
#         exhaustiveSearch(workdir, infodir, idents, outdir, params)
#     elif searchmode == SEARCH_MODE_DICT['kmercnt']:
#         quickSearch(workdir, infodir, idents, outdir)


# def exhaustiveSearch(workdir, infodir, idents, outdir, params):
#     for y in xrange(len(idents)):
#         for x in xrange(y+1, len(idents)):
#             pairoutDir = prepCSROutput(outdir, idents[y], idents[x])
#             mask1 = os.path.join(workdir, settings.MASKREPEATS_DIR, "%s_masked.fasta"%(idents[y]))
#             mask2 = os.path.join(workdir, settings.MASKREPEATS_DIR, "%s_masked.fasta"%(idents[x]))
#             subprocess.call(['Seanome.py',  '--threads', settings.WORK_THREADS, "seed_csr" , "-i1" , mask1, "-i2", mask2, "-o",  pairoutDir, "-l" , str(params.get('minlen',settings.MIN_SEQ_LEN)) , "-s" , str(params.get('sim', settings.MIN_SEQ_SIM))])
#     #logger.info(outdir)
#     #logger.info(os.path.join(infodir, settings.ORDER_FILE))
#     subprocess.call(['findBestOrder.py',  outdir, os.path.join(infodir, settings.ORDER_FILE)])


# def quickSearch(workdir, infodir, idents, outdir):
#     for name in idents:
#         outprefix = os.path.join(outdir, name)
#         subprocess.call(['findOrder.sh', os.path.join(workdir, settings.UCLUST_TO_SAM_DIR, '%s_pseudo_ref_parts.fasta'%(name)), outprefix, "2", "21"])
#     #logger.info(outdir)
#     #logger.info(os.path.join(infodir, settings.ORDER_FILE))
#     subprocess.call(['findBestOrder_quick.py',  outdir, '_kmer.filter', os.path.join(infodir, settings.ORDER_FILE)])


# class processCSR(threading.Thread):
#     def __init__(self, workdir, f, finalALNDir, MSARgnoutdir, MSARgncleanoutdir, MSARgntrimoutdir, SNPConsDir, SNPDir, SNPvcfDir, samSplitDir, semA, semB):
#         super(processCSR, self).__init__()
#         self.workdir = workdir
#         self.f = f
#         self.finalALNDir = finalALNDir
#         self.MSARgnoutdir = MSARgnoutdir
#         self.MSARgncleanoutdir = MSARgncleanoutdir 
#         self.MSARgntrimoutdir = MSARgntrimoutdir
#         self.SNPConsDir = SNPConsDir
#         self.SNPDir = SNPDir
#         self.SNPvcfDir = SNPvcfDir
#         self.semA = semA
#         self.semB = semB
#         self.samSplitDir = samSplitDir
#         self.start()

#     def run(self):
#         f = self.f
#         finalALNDir = self.finalALNDir
#         MSARgnoutdir = self.MSARgnoutdir
#         MSARgncleanoutdir = self.MSARgncleanoutdir
#         MSARgntrimoutdir = self.MSARgntrimoutdir
#         SNPConsDir = self.SNPConsDir
#         SNPDir = self.SNPDir
#         SNPvcfDir = self.SNPvcfDir
#         workdir = self.workdir
#         samSplitDir = self.samSplitDir
#         cleanid = os.path.splitext(f)[0]
#         convert = os.path.join(finalALNDir, f)
#         region = os.path.join(MSARgnoutdir, "%s.sam"%(f))
#         c_region = os.path.join(MSARgncleanoutdir, "%s.sam"%(f))
#         cb_region = os.path.join(MSARgncleanoutdir, "%s.bam"%(f))
#         trimout = os.path.join(MSARgntrimoutdir, cleanid)
#         consensus = os.path.join(SNPConsDir, f)
        
#         subprocess.call(['addConsensus.py', convert, convert])
#         with open('%s.log'%(trimout), "w") as logf:
#             subprocess.call(['trimal', '-in', convert, '-out', "%s.msa"%(trimout) , '-clustal', '-gt', '0.8', '-st', '0.001', '-cons', '60', '-colnumbering'], stdout = logf)
#         subprocess.call(['seqret', '-sequence', "%s.msa"%(trimout) , '-sformat1', 'clustal', '-outseq', "%s.fasta"%(trimout) , '-osformat2', 'fasta'])
#         subprocess.call(['samtools', 'faidx', "%s.fasta"%(trimout)])
#         subprocess.call(['Seanome.py',  '--threads', settings.WORK_THREADS , "inferSAM" , "-a" , convert, '-s', os.path.join(workdir, settings.UCLUST_TO_SAM_DIR), "-c", '-n', "-o", MSARgnoutdir])
#         subprocess.call(['cleanSAM.py',  '%s.msa'%(trimout,), region, "%s.log"%(trimout), c_region])
#         subprocess.call(['sam_to_bam.sh', c_region])
#         conseq = SeqIO.parse("%s.fasta"%(trimout), "fasta").next()
#         with open(consensus, "w") as o:
#             print >> o, ">%s\nN%s"%(conseq.id, str(conseq.seq))
#         subprocess.call(['samtools', 'faidx', consensus])
#         file_list = eval(subprocess.check_output(['split_sam.py', c_region, samSplitDir]) )
#         inputs = []
#         logger.info(file_list)
#         for ident, fpth in file_list:
#             bampth = "%s.bam"%(os.path.splitext(fpth)[0])
#             subprocess.call(['sam_to_bam.sh', fpth])
#             subprocess.call(['addReadGroup.sh', bampth, ident])
#             inputs.append(bampth)
#         subprocess.call(['runPlatypus.sh', ",".join(inputs), consensus, "%s.vcf"%(os.path.join(SNPvcfDir, f)), settings.WORK_THREADS])
#         subprocess.call(['vcfmod.py', '-f', '-b', cb_region  , '-i', '%s.vcf'%(os.path.join(SNPvcfDir, f)), '-o' , '%s.vcf'%(os.path.join(SNPDir, f)), '-p', '1', '-d', '%s.dat'%(os.path.join(SNPDir, f))])
#         #rdis_semaphore.release()    
#         self.semA.release()
#         self.semB.release()


# def executeMSABuilder(workdir, jid, idents, params):
#     SNPDir = os.path.join(workdir, settings.SNP_DIR)
#     createdir(SNPDir)
#     samSplitDir = os.path.join(SNPDir, "split_sam")
#     createdir(samSplitDir)
#     SNPConsDir = os.path.join(SNPDir, "consensus_fasta")
#     createdir(SNPConsDir)
#     SNPvcfDir = os.path.join(SNPDir, "called_vcf")
#     createdir(SNPvcfDir)

#     outdir = os.path.join(workdir, settings.MSA_DIR)
#     MSAInfoDir = os.path.join(outdir, 'information')
#     createdir(MSAInfoDir)
#     masked_files = dict([ (i, os.path.join(workdir, settings.MASKREPEATS_DIR, "%s_masked.fasta"%(i)), ) for i in idents])

#     generateOrdering(workdir, MSAInfoDir, idents, params, job.objects.get(pk = jid).searchMode)
    
#     ordering = eval(open(os.path.join(MSAInfoDir,  settings.ORDER_FILE)).read())
#     if len(ordering) < 2:
#         return
    
#     MSAoutdir = prepCSROutput(os.path.join(outdir, "01_seed_csr"), ordering[0], ordering[1])
#     subprocess.call(['Seanome.py',  '--threads', settings.WORK_THREADS, "seed_csr" , "-i1" , masked_files[ordering[0]], "-i2", masked_files[ordering[1]], "-o",  MSAoutdir, 
#                      "-l" , str(params.get('minlen',settings.MIN_SEQ_LEN)) , "-s" , str(params.get('sim', settings.MIN_SEQ_SIM))])
#     seedsdir = MSAoutdir

#     previous_full = MSAoutdir
#     MSAFoutdir = ""
#     previous = os.path.basename(MSAoutdir)
#     for i in ordering[2:]:
#         basedir = "%s_%s"%(previous, i)
#         MSAoutdir = os.path.join(outdir, "02_aln", basedir )
#         createdir(MSAoutdir)
#         MSAFoutdir = os.path.join(outdir, "02_fasta", basedir )
#         createdir(MSAFoutdir)
#         subprocess.call(['Seanome.py',  '--threads', settings.WORK_THREADS , "find_csr" , "-a" , previous_full, "-g", masked_files[i], "-o",  MSAoutdir, "-f", MSAFoutdir, 
#                          "-l" , str(params.get('minlen',settings.MIN_SEQ_LEN)) , "-s" , str(params.get('sim', settings.MIN_SEQ_SIM))])
#         previous = basedir
#         previous_full = MSAoutdir

#     finalFADir = MSAFoutdir
#     finalALNDir = previous_full

#     MSARgnoutdir = os.path.join(outdir, '04_infersam')
#     MSARgntrimoutdir = os.path.join(outdir, '03_trimAL')
#     MSARgncleanoutdir = os.path.join(outdir, '05_infersam_trimmed')
#     MSAconcatDir = os.path.join(outdir, '06_concat_trimeAL_msa')   
    
#     #pileoutdir = os.path.join(outdir, "07_trimmed_pileup" )
#     #createdir(pileoutdir)

#     createdir(MSARgnoutdir)
#     createdir(MSARgntrimoutdir)
#     createdir(MSARgncleanoutdir)
#     createdir(MSAconcatDir)
#     subprocess.call(['breakdown.py', finalFADir, os.path.join(MSAInfoDir, "find_csr") ] )
#     # DLS -- How can this be paralized to run X iterations of the loop at once?
#     # spin off X threads and wait for them to compelete..
#     threads = []
#     # used in an attempt to limit how many threads of processing a task spins up at any given time.
#     sem = threading.Semaphore(int(settings.REDIS_THREADS))
#     rdis_semaphore = Semaphore(Redis(), count = settings.REDIS_THREADS, namespace = 'SeanomeWorkers')
#     #sem = threading.Semaphore(6)
#     for f in os.listdir(finalALNDir):
#         rdis_semaphore.acquire()
#         sem.acquire()
#         threads.append(processCSR(workdir, f, finalALNDir, MSARgnoutdir, MSARgncleanoutdir, MSARgntrimoutdir, SNPConsDir, SNPDir, SNPvcfDir,  samSplitDir, rdis_semaphore, sem))
#     for t in threads:
#         t.join()
#     subprocess.call(['sambreakdown.py', MSARgnoutdir, os.path.join(MSAInfoDir, "pre_trim") ] )
#     subprocess.call(['sambreakdown.py', MSARgncleanoutdir, os.path.join(MSAInfoDir, "post_trim") ] )
    
#     for i in xrange(2, len(ordering)+1):
#         #logger.info(" ".join(["catALN.py", MSARgntrimoutdir, 'fasta', str(i), os.path.join(MSAInfoDir, settings.ORDER_FILE), os.path.join(MSAconcatDir, 'trimmed_cat_aln_%s'%(i))]))
#         subprocess.call(["catALN.py", MSARgntrimoutdir, 'fasta', str(i), os.path.join(MSAInfoDir, settings.ORDER_FILE), os.path.join(MSAconcatDir, 'trimmed_cat_aln_%s'%(i))])

#     return seedsdir, finalALNDir, finalFADir

# # --refFile == the MSA/fa with Con
# # --bamFiles == the sam (converted to bam) from the post trimAL 
# #python ~/programs/Platypus_0.5.2/Platypus.py callVariants --bamFiles=IH_20333_20467_WRG.bam,IWS_34172_34306_R_WRG.bam,Mincr_4999_5131_WRG.bam
# #--refFile=ref/outFile_9_146.fa --output=VariantCalls.vcf --filterDuplicates 0 --filterReadsWithUnmappedMates 0  --filterReadsWithDistantMates 0 --minMapQual 20   --minVarFreq 0.01 --minBaseQual 10 --genSNPs 1 --genIndels 0 -n 3 -p 2 --badReadsWindow 2  --maxSize 2 --assemblyRegionSize 2  --filterReadPairsWithSmallInserts 0 --maxVariants 100   --minGoodQualBases 10 --minFlank 2 --nCPU 2
