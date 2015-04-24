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
import yaml
import itertools
import cStringIO as StringIO

from Bio import AlignIO

# also in the utils folder but that is located somewhere else..
QUERY_CSR_AS_SEQS = """SELECT fileID, group_concat(seqID, '\t') AS IDs, group_concat(sequence, '\t') AS SEQS FROM csr GROUP BY fileID"""
QUERY_TRIMMED_CSR_AS_SEQS = """SELECT fileID, group_concat(seqID, '\t') AS IDs, group_concat(sequence, '\t') AS SEQS FROM trimmed_csr GROUP BY fileID"""
DOWNLOAD_TYPE = dict(post_trim_cat = 0, post_trim_aln = 1, post_trim_sam = 4, pre_trim_aln = 2, pre_trim_sam = 3, all = 5)

#http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks-in-python
def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def nparams(count):
    return ",".join("?"*count)

def getsqliteDB(toplvl):
    return sqlite3.connect(os.path.join(toplvl,  settings.SUFFIX_DB_LOCATION), check_same_thread = False)
 

def getSpecificFilesByMemberCount(con, count, markhasjson = False):
    if markhasjson:
        return con.execute("""SELECT B.name, A.fileID, count(*) as 'size', C.json FROM groups AS A JOIN files as B ON (A.fileID = B.id) LEFT JOIN trimmed_modvcf as C ON (C.fileID = B.id) GROUP BY A.fileID HAVING size = ?;""", (count,))
    else:
        return con.execute("""SELECT B.name, fileID, count(*) as 'size' FROM groups AS A JOIN files as B ON (A.fileID = B.id) GROUP BY fileID HAVING size = ?;""", (count,))


def getMemberCounts(con, sizeonly = False):
    if not sizeonly:
        return con.execute("""SELECT fileID, count(*) as 'size' FROM groups GROUP BY fileID;""")
    else:
        return con.execute("""SELECT count(*) as 'size' FROM groups GROUP BY fileID;""")

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
    username = "_".join([settings.SFTPUSER, str(uid), str(random.randint(0, 1000))])
    #username = "_".join([settings.SFTPUSER, str(uid)])
    cursor.execute("INSERT user (login_name, password) VALUES (%s, PASSWORD(%s));", [username , passwd])
    cursor.execute("""SELECT user_id FROM user where login_name = %s;""", [username])
    uid = cursor.fetchone()[0]    
    cursor.execute("UNLOCK TABLES;")
    return  username, passwd, uid
                        

def getUserName(uid):
    cursor = connections['sftponly'].cursor()
    cursor.execute("""SELECT login_name from user where user_id = %s;""", [uid])
    return cursor.fetchone()[0]
    

def removeSFTPUserbyName(username):
    cursor = connections['sftponly'].cursor()
    cursor.execute("""DELETE FROM user where login_name = %s;""", [username])


def removeSFTPUserbyUserID(uid):
    cursor = connections['sftponly'].cursor()
    cursor.execute("""DELETE FROM user where user_id = %s;""", [uid])


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
        return render_to_response(template, dict(jobj = jobj, min_sim = settings.MIN_SEQ_SIM, min_len = settings.MIN_SEQ_LEN ), context_instance =  RequestContext(request) )
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
    for i in xrange(len(flist)):
        if not params["N_%s"%(i)]:
            continue
        ident = "".join(params["N_%s"%(i)])
        ident = ident.lower().replace("_","-")
        idmap[ident].append("".join(params["F_%s"%(i)]))
        errors = [i for i,v in idmap.iteritems() if len(v) != 2]
    if not errors:
        jobj = job.objects.get(pk = jid)
        jobj.interrupt = interrupt
        jobj.stage = STAGES_REVERSE['queued']
        sparams = dict(samples = params.get("id_sample", 'multi'), asmprocess = params.get('id_asm', 'clust'), repeatmask = int(params.get('id_repeats', 1)), sim = float(params.get('id_sim',settings.MIN_SEQ_SIM)), minlen = int(params.get('id_slen', settings.MIN_SEQ_LEN)))
        jobj.scriptparams = json.dumps(sparams, separators=(',',':'))
        jobj.single = False if (params.get("id_sample", 'multi') == 'multi') else True
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
    config = yaml.load(open(os.path.join(jobfolder, "coverageLocation.yaml")))
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
    con = getsqliteDB(toplvl)
    graphdata = []
    for k, v in Counter([l[0] for l in getMemberCounts(con, True) ]).iteritems():
        graphdata.append( dict(label = k, v= int(v), value = max(math.log10(int(v)), math.log10(float(v)+0.5)), url = reverse("filterResults", args = [jid, k] ) ) )
    #if os.path.exists(os.path.join(toplvl, "csr", "search.order")):
    genomes = [ q.strip() for q in open(os.path.join(toplvl, "csr", "search.order")) ]
    numgenomes = len(genomes)
    if numgenomes == 1 and genomes[0] == 'combined':
        genomes_sorted = [("combined", sum([len(s[0])for s in con.execute("""SELECT sequence FROM trimmed_consensus;""") ]) )]
        numgenomes =  con.execute("""SELECT count(distinct species) FROM groups;""").fetchone()[0] 
    else:
        genomes_sorted = [ (i, len(SeqIO.parse(os.path.join(toplvl, i, "%s_pseudo_ref.fasta"%(i)), "fasta").next()) , ) for i in sorted(genomes)]

    return render_to_response(template, dict(jobj = jobj, numgenomes = numgenomes, genomelens = genomes_sorted, ordering = genomes, data = graphdata, tally= json.dumps(graphdata, separators=(',',':') ) ), context_instance = RequestContext(request) )


def ajaxGetSNPs(request, jid, ident):
    toplvl = os.path.join(settings.JOBDIR, jid)
    con = getsqliteDB(toplvl)
    jsondat = []
    for i in con.execute("""SELECT json FROM trimmed_modvcf WHERE fileID = ?;""", (ident,) ):
        jsondat = i[0]
    return HttpResponse(jsondat, content_type = "application/json")    

    
def ajaxGetMSA(request, jid, ident, clean):
    jobj = job.objects.get(pk = jid)
    toplvl = os.path.join(settings.JOBDIR, jid)
    con = getsqliteDB(toplvl)
    if jobj.single == False:
        qry =  QUERY_TRIMMED_CSR_AS_SEQS
        cons = "trimmed_consensus"
        if clean == '0':
            qry = QUERY_CSR_AS_SEQS
            cons = "consensus"
        data = con.execute("""SELECT A.fileID, B.sequence, A.IDs, A.SEQS FROM ( %(csr_as_seqeuences)s ) AS A JOIN %(contbl)s AS B ON (A.fileID = B.fileID) where A.fileID = ?;"""%dict(csr_as_seqeuences = qry, contbl = cons), (ident,) )
        rows = [  (r[0], ">%(seqID)s\n%(seq)s\n"%dict(seqID = "Consensus", seq = r[1]), "\n".join([">%s\n%s"%(i,s,)  for i, s in itertools.izip(r[2].split("\t") , r[3].split("\t") ) ] ), ) for r in data ]
        fastadata = "%s%s"%(rows[0][1], rows[0][2])
    else:
        data = con.execute("""SELECT B.fileID, B.sequence FROM trimmed_consensus AS B WHERE B.fileID = ?;""", (ident,) )
        rows = [  (r[0], ">%(seqID)s\n%(seq)s\n"%dict(seqID = "Consensus", seq = r[1]), ) for r in data ]
        fastadata = "%s"%(rows[0][1])        
    return HttpResponse(fastadata, content_type = "application/fasta")    


def viewAlignments(request, jid, count, template):
    jobj = job.objects.get(pk = jid)
    if jobj.stage != STAGES_REVERSE['done-s2']:
        return HttpResponseRedirect(reverse('runjob', args = [jid]))
    count = int(count)
    toplvl = os.path.join(settings.JOBDIR, jid)
    con = getsqliteDB(toplvl)
    # all the files that exist at a given species count...
    # 
    hasSNPs = []
    noSNPs = []
    for r in getSpecificFilesByMemberCount(con, count, True):
        if r[3]:
            hasSNPs.append( (r[0], r[1], True,) )
        else:
            noSNPs.append((r[0], r[1], False,) )
    
    tally = natsort.natsorted(hasSNPs) + natsort.natsorted(noSNPs) 
    #natsort.natsorted( [ (r[0], r[1], True if r[3] else False, ) for r in getSpecificFilesByMemberCount(con, count, True)] )
    #tally = natsort.natsorted( getSpecificFilesByMemberCount(con, count).fetchall() )
    return render_to_response(template, dict(jobj = jobj, count = count, showoptions = True, tally = tally), context_instance = RequestContext(request) )


def resultFilter(request, jid, count, template):
    jobj = job.objects.get(pk = jid)
    if jobj.stage != STAGES_REVERSE['done-s2']:
        return HttpResponseRedirect(reverse('runjob', args = [jid]))

    count = int(count)
    toplvl = os.path.join(settings.JOBDIR, jid)
    # get the file names of the csrs that have this particular
    con = getsqliteDB(toplvl)
    results = con.execute("""SELECT sum(coverage) as tot_cov,  count(*) as 'size' FROM groups GROUP BY fileID HAVING size = ?;""", (count,) )
    prebucket = Counter( (r[0] for r in results) )
    results = con.execute("""SELECT sum(trimmed_coverage) as tcov,  count(*) as 'size' FROM groups GROUP BY fileID HAVING size = ?;""", (count,) )
    postbucket = Counter( (r[0] for r in results) )

    total = 0
    prebrkcm = {}
    for k in sorted(prebucket.keys(), reverse = True):
        total += prebucket[k]
        prebrkcm[k] = total

    total = 0 
    postbrkcm = {}
    for k in sorted(postbucket.keys(), reverse = True):
        total += postbucket[k]
        postbrkcm[k] = total

    precovall = sorted(prebrkcm.iteritems(), key = lambda x: x[0])
    postcovall = sorted(postbrkcm.iteritems(), key = lambda x: x[0])

    return render_to_response(template, dict(jobj = jobj, count = count, showoptions = True, single = jobj.single,
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
    logfile = os.path.join(settings.JOBDIR, jid, "log.out")
    status = 0
    logdat = ""
    if not (jobj.stage == STAGES_REVERSE['queued'] or jobj.stage == STAGES_REVERSE['running']):
        status = 1
    else:
        logdat = "".join([i for i in open(logfile)][-50:])
        
    goto = reverse('runjob', args=[jid])
    return HttpResponse(json.dumps(dict(status=status, page=goto, logdat = logdat), separators=(',',':')), content_type = "application/json")


def downloadvcf(request, jid, count):
    jobj = job.objects.get(pk = jid)  
    count = int(count)
    toplvl = os.path.join(settings.JOBDIR, jid)
    con = getsqliteDB(toplvl)
    fids = tuple([str(f[1]) for f in getSpecificFilesByMemberCount(con, count)]) # fname and fid
    zstream = ZipFile(fileobj = None, compression = ZIP_DEFLATED)
    if jobj.single == False:
        for cfids in chunks(fids, 100):
            data = con.execute("""SELECT B.name, C.sequence, group_concat(A.seqID, '\t') AS IDs, group_concat(A.sequence, '\t') AS SEQS, D.vcf  
                          FROM trimmed_csr as A JOIN files as B ON (B.id = A.fileID) JOIN trimmed_consensus as C ON (C.fileID = A.fileID) JOIN trimmed_vcf AS D ON (D.fileID = A.fileID)
                          GROUP BY A.fileID HAVING A.fileID IN (""" + nparams(len(cfids)) + """);""", cfids )
            for f in data:
                fdat = "".join([ 
                        ">%(seqID)s\n%(seq)s\n"%dict(seqID = "Consensus", seq = f[1]),
                        "\n".join([">%s\n%s"%(i, s,)  for i, s in itertools.izip(f[2].split("\t") ,f[3].split("\t") ) ]) 
                        ])
                zstream.write(iterable = [fdat], arcname = os.path.join("post_trim_alignment", "%s.fasta"%(f[0]) ) )
                zstream.write(iterable = [f[-1]], arcname = os.path.join("SNPs", "%s.vcf"%(f[0])) )
    else:
        for cfids in chunks(fids, 100):
            data = con.execute("""SELECT B.name, C.sequence, D.vcf FROM files as B JOIN trimmed_consensus as C 
                          ON (C.fileID = B.id) JOIN trimmed_vcf AS D ON (D.fileID = B.id)
                          WHERE B.id IN (""" + nparams(len(cfids)) + """);""", cfids )
            for f in data:
                fdat =  ">%(seqID)s\n%(seq)s\n"%dict(seqID = "Consensus", seq = f[1])
                zstream.write(iterable = [fdat], arcname = os.path.join("post_trim_alignment", "%s.fasta"%(f[0]) ) )
                zstream.write(iterable = [f[-1]], arcname = os.path.join("SNPs", "%s.vcf"%(f[0])) )

    return streamArchive(zstream, "Job_%s_filter_%s_SNPs.zip"%(jid,count) )    
    


def MSADownload(request, jid, count, tid): #msa_download
    ### TODO: DLS - figure out a nice way to not iterator the sqlite execute command before dealing with the streaming

    jobj = job.objects.get(pk = jid)  
    tid = int(tid)
    count = int(count)

    toplvl = os.path.join(settings.JOBDIR, jid)
    zstream = ZipFile(fileobj = None, compression = ZIP_DEFLATED)

    if tid == DOWNLOAD_TYPE['post_trim_cat']: # the concat alns
        parentdir = os.path.join(toplvl, "concat_trimmed")
    
        zstream.write(os.path.join(parentdir,"%s_concat.msa"%(count) ), arcname = os.path.join("concat_aln", "%s_concat.msa"%(count)) )
        zstream.write(os.path.join(parentdir,"%s_concat.fasta"%(count)), arcname = os.path.join("concat_aln", "%s_concat.fasta"%(count)) )
        return streamArchive(zstream, "Job_%s_filter_%s_post-trim_concat.zip"%(jid, count) )

    # elif tid == DOWNLOAD_TYPE['all']:
    #     parentdir = os.path.join(toplvl, settings.MSA_DIR, "06_concat_trimeAL_msa")
    #     zstream.write(os.path.join(parentdir,"trimmed_cat_aln_%s.msa"%(count) ),os.path.join("post_trim_concat", "trimmed_cat_aln_%s.msa"%(count)) )
    #     zstream.write(os.path.join(parentdir,"trimmed_cat_aln_%s.fasta"%(count)),os.path.join("post_trim_concat", "trimmed_cat_aln_%s.fasta"%(count)) )
    #     parentaln = os.path.join(toplvl, settings.MSA_DIR, "03_trimAL")
    #     keep = [l.split()[0] for l in open(os.path.join(infodir, "find_csr.count")) if int(l.strip().split()[1]) == count]
    #     if os.path.exists(os.path.join(parentaln,"%s.msa.fasta"%(keep[0]))):
    #         for f in keep:
    #             zstream.write(os.path.join(parentaln,"%s.msa.fasta"%(f)),os.path.join("post_trim_alignment", "%s.fasta"%(f)) )
    #     else:
    #         for f in keep:
    #             zstream.write(os.path.join(parentaln,"%s.fasta"%(f)),os.path.join("post_trim_alignment", "%s.fasta"%(f)) )
    #     data = json.loads(jobj.data)
    #     parentaln = data['aln']
    #     for f in keep:
    #         zstream.write(os.path.join(parentaln,f),os.path.join("pre_trim_alignment", "%s.fasta"%(f)) )
    #     parentaln = os.path.join(toplvl, settings.MSA_DIR, "04_infersam")
    #     for f in keep:
    #         zstream.write(os.path.join(parentaln,"%s.sam"%(f)),os.path.join("pre_trim_alignment_sam", "%s.sam"%(f)) )
    #     parentaln = os.path.join(toplvl, settings.MSA_DIR, "05_infersam_trimmed")
    #     for f in keep:
    #         zstream.write(os.path.join(parentaln,"%s.sam"%(f)),os.path.join("post_trim_alignment_sam", "%s.sam"%(f)) )
    #     parentsnp = os.path.join(toplvl, "15_SNP_calling")
    #     for f in keep:
    #         zstream.write(os.path.join(parentsnp,"%s.vcf"%(f)),os.path.join("post_trim_SNPs", "%s.vcf"%(f)) )
    #     return streamArchive(zstream, "Job_%s_filter_%s_all.zip"%(jid,count) )

   
    mincov = -1
    if request.method == 'GET':
        mincov = int(request.GET.get('cutoff', -1))
    con = getsqliteDB(toplvl)

    results_pretrim = tuple([ str(r[0]) for r in con.execute("""SELECT fileID, sum(coverage) as tcov, count(*) as 'size' FROM groups GROUP BY fileID HAVING size = ? AND tcov >= ?;""", (count, mincov,)) ])
    results_posttrim = tuple([ str(r[0]) for r in con.execute("""SELECT fileID, sum(trimmed_coverage) as tcov, count(*) as 'size' FROM groups GROUP BY fileID HAVING size = ? AND tcov >= ?;""", (count, mincov,)) ])
    
    if tid == DOWNLOAD_TYPE['post_trim_aln']: # non cat alns in zip, post-trim
        if jobj.single == False:
            for cdata in chunks(results_posttrim, 100):
                data = con.execute("""SELECT B.name, C.sequence, group_concat(A.seqID, '\t') AS IDs, group_concat(A.sequence, '\t') AS SEQS 
                              FROM trimmed_csr as A JOIN files as B ON (B.id = A.fileID) JOIN trimmed_consensus as C ON (C.fileID = A.fileID) 
                              GROUP BY A.fileID HAVING A.fileID IN (""" + nparams(len(cdata)) + """);""", cdata )
                for f in data:
                    fdat = "".join([ 
                            ">%(seqID)s\n%(seq)s\n"%dict(seqID = "Consensus", seq = f[1]),
                            "\n".join([">%s\n%s"%(i, s,)  for i, s in itertools.izip(f[2].split("\t"), f[3].split("\t") ) ]) 
                            ])
                    zstream.write(iterable = [fdat], arcname = os.path.join("post_trim_alignment", "%s.fasta"%(f[0]) ) )
        else:
            for cdata in chunks(results_posttrim, 100):
                data = con.execute("""SELECT B.name, C.sequence
                              FROM files as B JOIN trimmed_consensus as C ON (C.fileID = B.id) 
                              WHERE B.id IN (""" + nparams(len(cdata)) + """);""", cdata )
                for f in data:
                    fdat = ">%(seqID)s\n%(seq)s\n"%dict(seqID = "Consensus", seq = f[1])
                    zstream.write(iterable = [fdat], arcname = os.path.join("post_trim_alignment", "%s.fasta"%(f[0]) ) )
        return streamArchive(zstream, "Job_%s_filter_%s_post-trim_aln_mincov_%s.zip"%(jid, count, mincov) )

    elif tid == DOWNLOAD_TYPE['post_trim_sam']: # non cat alns in zip. pre-trim        
        for cdata in chunks(results_posttrim, 100):
            data = con.execute("""SELECT A.name, B.sam FROM files as A JOIN trimmed_inferSAM as B ON (A.id = B.fileID) WHERE A.id IN (""" + nparams(len(cdata)) + """);""", cdata)
            for f in data:
                zstream.write(iterable = [f[1]], arcname = os.path.join("post_trim_alignment_sam", "%s.sam"%(f[0]) ) )
        return streamArchive(zstream, "Job_%s_filter_%s_pre-trim_aln_mincov_%s.zip"%(jid, count, mincov) )

    elif tid == DOWNLOAD_TYPE['pre_trim_aln']: #pre-trim sam
        if jobj.single == False:
            for cdata in chunks(results_pretrim, 100):
                data = con.execute("""SELECT B.name, C.sequence, group_concat(A.seqID, '\t') AS IDs, group_concat(A.sequence, '\t') AS SEQS 
                              FROM csr AS A JOIN files AS B ON (B.id = A.fileID) JOIN consensus as C ON (C.fileID = A.fileID) 
                              GROUP BY A.fileID HAVING A.fileID IN (""" + nparams(len(cdata)) + """);""",  cdata )
                for f in data:
                    fdat = "".join([ 
                            ">%(seqID)s\n%(seq)s\n"%dict(seqID = "Consensus", seq = f[1]),
                            "\n".join([">%s\n%s"%(i, s,)  for i, s in itertools.izip(f[2].split("\t") ,f[3].split("\t") ) ]) 
                            ])
                    zstream.write(iterable = [fdat], arcname = os.path.join("pre_trim_alignment", "%s.fasta"%(f[0]) ) )
        else:
            for cdata in chunks(results_pretrim, 100):
                data = con.execute("""SELECT B.name, C.sequence
                              FROM files as B JOIN trimmed_consensus as C ON (C.fileID = B.id) 
                              WHERE B.id IN (""" + nparams(len(cdata)) + """);""",  cdata )
                for f in data:
                    fdat =  ">%(seqID)s\n%(seq)s\n"%dict(seqID = "Consensus", seq = f[1])
                    zstream.write(iterable = [fdat], arcname = os.path.join("pre_trim_alignment", "%s.fasta"%(f[0]) ) )
        return streamArchive(zstream, "Job_%s_filter_%s_pre-trim_sam_mincov_%s.zip"%(jid, count, mincov) )

    elif tid == DOWNLOAD_TYPE['pre_trim_sam']: #post-trim sam
        for cdata in chunks(results_pretrim, 100):
            data = con.execute("""SELECT A.name, B.sam FROM files as A JOIN inferSAM as B ON (A.id = B.fileID) WHERE A.id IN (""" + nparams(len(cdata)) + """);""", cdata)
            for f in data:
                zstream.write(iterable = [f[1]] ,arcname = os.path.join("pre_trim_alignment_sam", "%s.sam"%(f[0]) ) )        
        return streamArchive(zstream, "Job_%s_filter_%s_post-trim_sam_mincov_%s.zip"%(jid, count, mincov) )
