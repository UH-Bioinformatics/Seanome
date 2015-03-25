from django.template import RequestContext
from django.http import HttpResponseRedirect, HttpResponse, StreamingHttpResponse

from django.shortcuts import render_to_response, render

from django.conf import settings
from django.core.urlresolvers import reverse
from django.core.exceptions import ObjectDoesNotExist
from django.utils.encoding import smart_str

from django.db.models import Q, Max, Count
from django.db import connections

from django.core.cache import cache
from django.views.decorators.cache import cache_page, cache_control

from .filetransfer import *
from .models import *
from .forms import *
from .tasks import *
import natsort
import json
import tempfile
import os
import random
import string
import subprocess
import math
from collections import defaultdict, Counter
from zipstream import *
import sqlite3 


# Create your views here.
def createSFTPUser():
    # TODO: need to figure out a neat, clean way to generate a user name that wont have collisions!
    # DLS: attempted to work around this by using a table lock 
    length = 16
    chars = string.ascii_letters + string.digits + '!@#$%^&*()-_=+/\{}[]'   
    passwd = ''.join(random.choice(chars) for i in xrange(length))
    cursor = connections['sftponly'].cursor()
    cursor.execute("""LOCK TABLE user WRITE, user as userRead READ;""")
    cursor.execute("""SELECT COUNT(*) FROM user;""")
    uid = cursor.fetchone()[0]
    username = "_".join([settings.SFTPUSER, str(uid) ,str(random.randint(0, 1000))])
    #username = "_".join([settings.SFTPUSER, str(uid)])
    cursor.execute("INSERT user (login_name, password) VALUES (%s, PASSWORD(%s));", [username , passwd])
    cursor.execute("""SELECT user_id FROM user where login_name = %s;""", [username])
    uid = cursor.fetchone()[0]    
    cursor.execute("UNLOCK TABLES;")
    return  username, passwd, uid
                        

def getUserName(uid):
    cursor = connections['sftponly'].cursor()
    cursor.execute("""SELECT login_name from user where user_id = %s;""",[uid])
    return cursor.fetchone()[0]
    

def removeSFTPUserbyName(username):
    cursor = connections['sftponly'].cursor()
    cursor.execute("""DELETE FROM user where login_name = %s;""",[username])


def removeSFTPUserbyUserID(uid):
    cursor = connections['sftponly'].cursor()
    cursor.execute("""DELETE FROM user where user_id = %s;""",[uid])


def beginJob(request, template):
    if request.method == 'POST':
        form = JobStarter(request.POST)
        if form.is_valid():
            username, passwd, uid = createSFTPUser()
            form = None
            address = settings.SERVERADDRESS            
            jobsubmission = job(sftp_id = uid)
            jobsubmission.save()
            return render_to_response(template, dict(usr = username, passwd = passwd, address = address, job = jobsubmission.id), context_instance =  RequestContext(request) )
    else:
        form = JobStarter()
    return render_to_response(template, dict(form = form), context_instance =  RequestContext(request) )


def ajaxRefreshFileDir(request, jid):
    if not request.is_ajax():
        return
    jobj = job.objects.get(pk = jid)
    userhome = getUserName(jobj.sftp_id)
    lstout = listdir(userhome)
    return HttpResponse(json.dumps(lstout, separators=(',',':')), content_type = "application/json")

    
def listdir(userhome):
    orighome = os.path.join(settings.UPLOADED, userhome)
    try:
        lstout = subprocess.check_output(['sudo', '/bin/ls', '-1',orighome])
        lstout = [ l for l in lstout.strip().split("\n") if l.strip()]
    except:
        lstout = []
    return lstout


def prepstage(request, jid, jobj, template):
    if request.method != 'POST':        
        return render_to_response(template, dict(jobj = jobj), context_instance =  RequestContext(request) )
    userhome = getUserName(jobj.sftp_id)
    orighome = os.path.join(settings.UPLOADED, userhome)
    newhome = os.path.join(settings.JOBDIR, jid)
    flist = listdir(userhome)
    os.system("""sudo /bin/chown -R celery:www-data %s"""%(orighome))
    os.system("""sudo /bin/mv %s %s"""%(orighome, newhome))       
    removeSFTPUserbyUserID(jobj.sftp_id)
    params = request.POST
    idmap = defaultdict(list)
    errors = None
    interrupt = params.get('id_interrupt', "Continue")
    if interrupt == 'Interrupt':
        interrupt = True
    else:
        interrupt = False
    
    per = 2
    #if params.get('id_asm', 'clust') == 'preasm':
    #    per = 3

    for i in xrange(len(flist)):
        if not params["N_%s"%(i)]:
            continue
        ident = "".join(params["N_%s"%(i)])
        ident = ident.lower().replace("_","-")
        idmap[ident].append("".join(params["F_%s"%(i)]))
        errors = [i for i,v in idmap.iteritems() if len(v) != per]
    if not errors:
        jobj = job.objects.get(pk = jid)
        jobj.interrupt = interrupt
        jobj.stage = STAGES_REVERSE['queued']
        sparams = dict(asmprocess = params.get('id_asm', 'clust'), repeatmask = int(params.get('id_repeats', 1)), sim = float(params.get('id_sim', 0.94)), minlen = int(params.get('id_slen', 150)))
        jobj.scriptparams = json.dumps(sparams, separators=(',',':'))
        jobj.save()
        seanomeGO.delay(newhome, jid, idmap)
    else:
        jobj = job.objects.get(pk = jid)
        jobj.stage = STAGES_REVERSE['canceled']
        jobj.save()
        seanomeBad.delay(newhome)
    return render_to_response("jobWaiting.html", dict(jobj = jobj, errors = errors), context_instance = RequestContext(request) )   


def filterLevel(request, jid, jobj, template):
    jobfolder = os.path.join(settings.JOBDIR, jid)
    if request.method == 'POST':
        ranged = defaultdict(list)
        for k,v in request.POST.iteritems():
            if not k.startswith("data_"):
                continue
            k = "_".join(k.split("_")[1:-1])
            ranged[k].append(int(v))
        jobj = job.objects.get(pk = jid)
        jobj.snpcaller = request.POST.get('id_snper', 1) 
        jobj.searchMode =  request.POST.get('id_ordering',2)
        jobj.conMode =  request.POST.get('id_consen', 1)
        jobj.stage = STAGES_REVERSE['queued']
        jobj.save()
        seanomeGOstage2.delay(jobfolder, jid, ranged)
        return None 
    filestats = []
    maxsize = 0
    config = yaml.load(open(os.path.join(workdir, "config.yaml")))
    for f in config['coverage']:
        jstr = json.loads(open(f).read())
        convert = []
        total = 0
        for i in jstr[1][::-1]:
            total += i[1]
            convert.append([i[0],total, i[1]])            
        jstr[1] = convert[::-1]
        filestats.append(jstr)
    jstr = json.dumps(filestats, separators=(',',':'))
    return render_to_response(template, dict(jobj = jobj, data = jstr), context_instance = RequestContext(request) )   


def runJob(request, jid, template):
    jobj = job.objects.get(pk = jid)
    if jobj.stage == STAGES_REVERSE['prep']:
        return prepstage(request, jid, jobj, "jobRunner.html")
    elif jobj.stage == STAGES_REVERSE['queued'] or jobj.stage == STAGES_REVERSE['running']:
        return render_to_response( "jobWaiting.html", dict(jobj = jobj), context_instance = RequestContext(request) )   
    elif jobj.stage == STAGES_REVERSE['done-s1']:
            result = filterLevel(request, jid, jobj, "filterLevel.html") 
            if result:
                return result 
            else:
                jobj = job.objects.get(pk = jid)
                template = "jobWaiting.html"
    elif jobj.stage == STAGES_REVERSE['done-s2']:
        # probably redirect to another URL?
        return HttpResponseRedirect(reverse('completedjob', args = [jid]))
    elif jobj.stage == STAGES_REVERSE['canceled']:
        return render_to_response( "jobWaiting.html", dict(jobj = jobj), context_instance = RequestContext(request) )   
    else:
        # catch all of other invalid stage values
        pass
    return render_to_response(template, dict(jobj = jobj), context_instance = RequestContext(request) )   


def completedJob(request, jid, template):
    jobj = job.objects.get(pk = jid)
    if jobj.stage != STAGES_REVERSE['done-s2']:
        return HttpResponseRedirect(reverse('runjob', args = [jid]))

    toplvl = os.path.join(settings.JOBDIR, jid)
    con = sqlite3.connect(os.path.join(toplvl, "seanome.db3"), check_same_thread=False)
    
    infodir = os.path.join(toplvl, settings.MSA_DIR, "information")
    tally = []
    for k, v in Counter([l[1] for l in con.execute("""SELECT fileID, count(*) as 'size' FROM trimmed_csr GROUP BY fileID;""")]).iteritems():
        tally.append( dict(label = k, v= int(v), value = max(math.log10(int(v)), math.log10(float(v)+0.5)), url = reverse("filterResults", args = [jid, k] ) ) )
    genomes = [ q.strip() for q in open(os.path.join(toplvl, "csr", "serach.order")) ]
    numgenomes = len(genomes)
    genomes_sorted = [ (i, len(SeqIO.parse(os.path.join(toplvl, settings.UCLUST_TO_SAM_DIR, "%s_pseudo_ref.fasta"%(i)), "fasta").next()) , ) for i in sorted(genomes)]
    # we have finished everythign.. we need to collect some data to display out to the user...
    return render_to_response(template, dict(jobj = jobj, numgenomes = numgenomes, genomelens = genomes_sorted, ordering = genomes, data = tally, tally= json.dumps(tally, separators=(',',':') ) ), context_instance = RequestContext(request) )


def ajaxGetSNPs(request, jid, ident):
    toplvl = os.path.join(settings.JOBDIR, jid)
    SNPdir = os.path.join(toplvl, settings.SNP_DIR, "%s.dat"%(ident))
    return HttpResponse(open(SNPdir).read(), content_type = "application/json")    

    
def ajaxGetMSA(request, jid, ident, clean):
    #if not request.is_ajax():
    #    return
    toplvl = os.path.join(settings.JOBDIR, jid)
    MSAdir = os.path.join(toplvl, settings.MSA_DIR)
    if clean == '1':
        #MSAdir = os.path.join(MSAdir, "07_trimmed_pileup", "%s.fasta"%(ident))
        #if not os.path.exists(MSAdir):
        MSAdir = os.path.join(MSAdir, "03_trimAL", "%s.fasta"%(ident))
        if not os.path.exists(MSAdir):
            MSAdir = os.path.join(os.path.join(toplvl, settings.MSA_DIR), "03_trimAL", "%s.msa.fasta"%(ident))
        
    elif clean == '0':
        jobj = job.objects.get(pk = jid)
        data = json.loads(jobj.data)
        MSAdir = os.path.join(data['aln'], ident)
    return HttpResponse(open(MSAdir).read(), content_type = "application/fasta")    



def viewAlignments(request, jid, count, template):
    jobj = job.objects.get(pk = jid)
    if jobj.stage != STAGES_REVERSE['done-s2']:
        return HttpResponseRedirect(reverse('runjob', args = [jid]))

    count = int(count)
    toplvl = os.path.join(settings.JOBDIR, jid)
    infodir = os.path.join(toplvl, settings.MSA_DIR, "information")
    tally = [ l.split()[0] for l in open(os.path.join(infodir, "find_csr.count")) if int(l.strip().split()[1]) == count]
    tally = natsort.natsorted(tally)
    return render_to_response(template, dict(jobj = jobj, count = count, showoptions = True,
                                             tally = tally), context_instance = RequestContext(request) )


def resultFilter(request, jid, count, template):
    jobj = job.objects.get(pk = jid)
    if jobj.stage != STAGES_REVERSE['done-s2']:
        return HttpResponseRedirect(reverse('runjob', args = [jid]))

    count = int(count)
    toplvl = os.path.join(settings.JOBDIR, jid)
    infodir = os.path.join(toplvl, settings.MSA_DIR, "information")

    tally = [ l.split()[0] for l in open(os.path.join(infodir, "find_csr.count")) if int(l.strip().split()[1]) == count]

    stally = set(tally)
    prebucket = defaultdict(int)
    postbucket = defaultdict(int)
    for pre, post in zip(open(os.path.join(infodir, "pre_trim.cov")), open(os.path.join(infodir, "post_trim.cov")) ):
        pre = pre.split("\t")
        post = post.split("\t")
        if pre[0] in stally:
            v = int(json.loads(pre[1])[0][1])
            prebucket[v] += 1
        if post[0] in stally:
            v = int(json.loads(post[1])[0][1])
            postbucket[v] += 1

    total = 0
    prebrkcm = {}
    for k in sorted(prebucket.keys(), reverse = True):
        total += prebucket[k]
        prebrkcm[k] = total

    total =0 
    postbrkcm = {}
    for k in sorted(postbucket.keys(), reverse = True):
        total += postbucket[k]
        postbrkcm[k] = total
    o = open("/tmp/data", "w")
    o.write(str(postbrkcm))
    o.close()

    precovall = sorted(prebrkcm.iteritems(), key = lambda x: x[0])
    postcovall = sorted(postbrkcm.iteritems(), key = lambda x: x[0])

    #precovall = sorted(prebucket.iteritems(), key = lambda x: x[0])
    #postcovall = sorted(postbucket.iteritems(), key = lambda x: x[0])

    return render_to_response(template, dict(jobj = jobj, count = count, showoptions = True,
                                             #tally = tally,
                                             prebreak = json.dumps( precovall, separators=(',',':') ),
                                             postbreak = json.dumps( postcovall, separators=(',',':') ),
                                             prebrkcm = json.dumps( prebrkcm, separators=(',',':') ),
                                             postbrkcm = json.dumps(postbrkcm, separators=(',',':') ),
                                             precovall = precovall, postcovall = postcovall ), context_instance = RequestContext(request) )


def cancelJob(request, jid, template):
    jobj = job.objects.get(pk = jid)
    if jobj.stage != 0:
         return HttpResponseRedirect(reverse('startjob'))
    if jobj.stage == 0 and request.method == 'POST':
        try:
            userhome = getUserName(jobj.sftp_id)
            orighome = os.path.join(settings.UPLOADED, userhome)
            os.system("""sudo /bin/chown -R celery:www-data %s"""%(orighome))
            removeSFTPUserbyUserID(jobj.sftp_id)
            delhome = os.path.join(settings.UPLOADED, jid)
            os.system("""sudo /bin/mv %s %s"""%(orighome, delhome))
            os.system("""sudo /bin/rm -rf %s"""%(delhome))
        except:
            pass
        jobj.stage = STAGES_REVERSE['canceled'] 
        jobj.save()
        return HttpResponseRedirect(reverse('startjob'))
    return render_to_response(template, dict(jobj = jobj), context_instance = RequestContext(request) )        

        
def staticPage(request, template):
    return render_to_response(template, dict(msg = "Hello World!"), context_instance = RequestContext(request) )


def ajaxJobStatus(request, jid):
    if not request.is_ajax():
        return
    jobj = job.objects.get(pk = jid)
    status = 0
    if not (jobj.stage == STAGES_REVERSE['queued'] or jobj.stage == STAGES_REVERSE['running']):
        status = 1
    goto = reverse('runjob', args=[jid])
    return HttpResponse(json.dumps(dict(status=status, page=goto), separators=(',',':')), content_type = "application/json")


def downloadvcf(request, jid, count):
    jobj = job.objects.get(pk = jid)  
    count = int(count)
    toplvl = os.path.join(settings.JOBDIR, jid)
    infodir = os.path.join(toplvl, settings.MSA_DIR, "information")
    zstream = ZipFile(fileobj = None, compression = ZIP_DEFLATED)
    tally = [ l.split()[0] for l in open(os.path.join(infodir, "find_csr.count")) if int(l.strip().split()[1]) == count]
    parentsnp = os.path.join(toplvl, "15_SNP_calling")
    parentaln = os.path.join(toplvl,  settings.MSA_DIR,  "03_trimAL")
    for f in tally:
        zstream.write(os.path.join(parentsnp,"%s.vcf"%(f)),os.path.join("SNPs", "%s.vcf"%(f)) )
        zstream.write(os.path.join(parentaln,"%s.fasta"%(f)),os.path.join("post_trim_alignment", "%s.fasta"%(f)) )
        #zstream.write(os.path.join(parentaln,"%s.msa"%(f)),os.path.join("aln_msa", "%s.msa"%(f)) )
    return streamArchive(zstream, "Job_%s_filter_%s_SNPs.zip"%(jid,count) )    
    

DOWNLOAD_TYPE=dict(post_trim_cat=0, post_trim_aln=1, post_trim_sam=2, pre_trim_aln=3, pre_trim_sam=4, all=5)
def MSADownload(request, jid, count, tid): #msa_download
    jobj = job.objects.get(pk = jid)  
    tid = int(tid)
    count = int(count)

    toplvl = os.path.join(settings.JOBDIR, jid)
    infodir = os.path.join(toplvl, settings.MSA_DIR, "information")
    zstream = ZipFile(fileobj = None, compression = ZIP_DEFLATED)
    if tid == DOWNLOAD_TYPE['post_trim_cat']: # the concat alns
        parentdir = os.path.join(toplvl, settings.MSA_DIR, "06_concat_trimeAL_msa")
        zstream.write(os.path.join(parentdir,"trimmed_cat_aln_%s.msa"%(count) ),os.path.join("post_trim_concat", "trimmed_cat_aln_%s.msa"%(count)) )
        zstream.write(os.path.join(parentdir,"trimmed_cat_aln_%s.fasta"%(count)),os.path.join("post_trim_concat", "trimmed_cat_aln_%s.fasta"%(count)) )
        return streamArchive(zstream, "Job_%s_filter_%s_post-trim_concat.zip"%(jid,count) )
    elif tid == DOWNLOAD_TYPE['all']:
        parentdir = os.path.join(toplvl, settings.MSA_DIR, "06_concat_trimeAL_msa")
        zstream.write(os.path.join(parentdir,"trimmed_cat_aln_%s.msa"%(count) ),os.path.join("post_trim_concat", "trimmed_cat_aln_%s.msa"%(count)) )
        zstream.write(os.path.join(parentdir,"trimmed_cat_aln_%s.fasta"%(count)),os.path.join("post_trim_concat", "trimmed_cat_aln_%s.fasta"%(count)) )
        parentaln = os.path.join(toplvl, settings.MSA_DIR, "03_trimAL")
        keep = [l.split()[0] for l in open(os.path.join(infodir, "find_csr.count")) if int(l.strip().split()[1]) == count]
        if os.path.exists(os.path.join(parentaln,"%s.msa.fasta"%(keep[0]))):
            for f in keep:
                zstream.write(os.path.join(parentaln,"%s.msa.fasta"%(f)),os.path.join("post_trim_alignment", "%s.fasta"%(f)) )
        else:
            for f in keep:
                zstream.write(os.path.join(parentaln,"%s.fasta"%(f)),os.path.join("post_trim_alignment", "%s.fasta"%(f)) )
        data = json.loads(jobj.data)
        parentaln = data['aln']
        for f in keep:
            zstream.write(os.path.join(parentaln,f),os.path.join("pre_trim_alignment", "%s.fasta"%(f)) )
        parentaln = os.path.join(toplvl, settings.MSA_DIR, "04_infersam")
        for f in keep:
            zstream.write(os.path.join(parentaln,"%s.sam"%(f)),os.path.join("pre_trim_alignment_sam", "%s.sam"%(f)) )
        parentaln = os.path.join(toplvl, settings.MSA_DIR, "05_infersam_trimmed")
        for f in keep:
            zstream.write(os.path.join(parentaln,"%s.sam"%(f)),os.path.join("post_trim_alignment_sam", "%s.sam"%(f)) )
        parentsnp = os.path.join(toplvl, "15_SNP_calling")
        for f in keep:
            zstream.write(os.path.join(parentsnp,"%s.vcf"%(f)),os.path.join("post_trim_SNPs", "%s.vcf"%(f)) )
        return streamArchive(zstream, "Job_%s_filter_%s_all.zip"%(jid,count) )


    
    mincov = -1
    if request.method == 'GET':
        mincov = int(request.GET.get('cutoff', -1))

    tally = set( ( l.split()[0] for l in open(os.path.join(infodir, "find_csr.count")) if int(l.strip().split()[1]) == count) )
    flist = []
    fcov = "post_trim.cov"
    if tid == DOWNLOAD_TYPE['pre_trim_sam'] or tid == DOWNLOAD_TYPE['pre_trim_aln']:
        fcov = "pre_trim.cov"
   
    for l in open(os.path.join(infodir, fcov)):
        f, info = l.split("\t")
        if f not in tally:
            continue
        cov = int(json.loads(info)[0][1])
        if cov >= mincov:
            flist.append(f)

    if tid == DOWNLOAD_TYPE['post_trim_aln']: # non cat alns in zip, post-trim
        parentaln = os.path.join(toplvl, settings.MSA_DIR, "03_trimAL")
        if os.path.exists(os.path.join(parentaln,"%s.msa.fasta"%(flist[0]))):
            for f in flist:
                zstream.write(os.path.join(parentaln,"%s.msa.fasta"%(f)),os.path.join("post_trim_alignment", "%s.fasta"%(f)) )
        else:
            for f in flist:
                zstream.write(os.path.join(parentaln,"%s.fasta"%(f)),os.path.join("post_trim_alignment", "%s.fasta"%(f)) )
        return streamArchive(zstream, "Job_%s_filter_%s_post-trim_aln_mincov_%s.zip"%(jid,count, mincov) )
    elif tid == DOWNLOAD_TYPE['post_trim_sam']: # non cat alns in zip. pre-trim
        data = json.loads(jobj.data)
        parentaln = data['aln']
        for f in flist:
            zstream.write(os.path.join(parentaln,f),os.path.join("pre_trim_alignment", "%s.fasta"%(f)) )
        return streamArchive(zstream, "Job_%s_filter_%s_pre-trim_aln_mincov_%s.zip"%(jid,count,mincov) )
    elif tid == DOWNLOAD_TYPE['pre_trim_aln']: #pre-trim sam
        parentaln = os.path.join(toplvl, settings.MSA_DIR, "04_infersam")
        for f in flist:
            zstream.write(os.path.join(parentaln,"%s.sam"%(f)),os.path.join("pre_trim_alignment_sam", "%s.sam"%(f)) )
        return streamArchive(zstream, "Job_%s_filter_%s_pre-trim_sam_mincov_%s.zip"%(jid,count,mincov) )
    elif tid == DOWNLOAD_TYPE['pre_trim_sam']: #post-trim sam
        parentaln = os.path.join(toplvl, settings.MSA_DIR, "05_infersam_trimmed")
        for f in flist:
            zstream.write(os.path.join(parentaln,"%s.sam"%(f)),os.path.join("post_trim_alignment_sam", "%s.sam"%(f)) )
        return streamArchive(zstream, "Job_%s_filter_%s_post-trim_sam_mincov_%s.zip"%(jid,count,mincov) )
