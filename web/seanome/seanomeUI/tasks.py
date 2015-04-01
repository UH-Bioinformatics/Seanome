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


@shared_task
def seanomeGO(workdir, jid, idmap):
    jobj = job.objects.get(pk = jid)

    jobj.stage = STAGES_REVERSE['running']
    sparams = json.loads(jobj.scriptparams)
    jobj.save()

    os.system("""chmod -R 755 %s"""%(workdir))
    setupUserEnv(workdir)
    close_old_connections()

    # just in case something very odd occurs and the script fails to make the directories...
    os.system("""chmod -R 755 %s"""%(workdir))

    
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
       

   if(sparams.get('repeatmask', 0) == 1):
       os.system("SeanomeWrapper.py -c %s -d %s -t %s -j 1 multiple"%( os.path.join(workdir, "config.yaml") , os.path.join(workdir, "seanome.db3") , settings.WORK_THREADS))
   else:
       os.system("SeanomeWrapper.py -c %s -d %s -t %s -j 1 -s multiple"%( os.path.join(workdir, "config.yaml") , os.path.join(workdir, "seanome.db3") , settings.WORK_THREADS))

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
    jobj = job.objects.get(pk = jid)
    sparams = json.loads(jobj.scriptparams)
    jobj.stage = STAGES_REVERSE['running']
    jobj.save()
    setupUserEnv(workdir)
    close_old_connections()

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
                os.system("""%s -c %s -z %s  >> log.out 2>> log.err """%(l, v[0], v[1]))
            else:
                os.system("""%s %s %s  >> log.out 2>> log.err """%(l, 3, 200))
        rdis_semaphore.release()

    rdis_semaphore.acquire()
    os.system("""bash child_runner_script_3.sh -l %s -s %s >> log.out 2>> log.err """%( str(params.get('minlen',settings.MIN_SEQ_LEN)) , str(params.get('sim', settings.MIN_SEQ_SIM) )
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

