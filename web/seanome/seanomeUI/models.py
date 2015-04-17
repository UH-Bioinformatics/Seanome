from django.db import models
from django.contrib import admin

# Create your models here.

STAGES=((-1, 'canceled',), (0, 'prep',), (1, 'queued',), (2, 'running',),  (3, 'done-s1',), (4, 'done-s2',))
STAGES_DICT = dict(STAGES)
STAGES_REVERSE= dict([ (v[1],v[0],) for v in STAGES])

SEARCH_MODES = ((2, 'Seed and order (long run time)',), (1, 'Kmer counting (quick)') ) 
SEARCH_MODE_DICT = dict(exhaustive = 2, kmercnt = 1)

SNP_MODES = ((1, 'Platypus',), (2, 'Free bayes'), (3, 'GATK'), ) 
SNP_MODE_DICT = dict([ (v[1], v[0],) for v in SNP_MODES])

CON_MODE = ((1, 'IUPAC using majority',), (2, 'other'))
CON_MODE_DICT = dict(IUPAC_M = 1, OTHER=2)

BUILD_MODE = ((1,'Clustering',),(2, 'Assemble',),)
FILTER_REPEATS = ((0, 'No',),(1, 'Yes',),)

class job(models.Model):
    id = models.AutoField(primary_key = True)
    sftp_id = models.IntegerField(blank = False, null = False, db_index = True)
    stage = models.IntegerField(blank = False, null = False, choices = STAGES, default = 0)
    data = models.TextField(blank = True, null = True)
    interrupt = models.BooleanField(default = True)
    searchMode = models.IntegerField(default = 1, choices = SEARCH_MODES)
    snpcaller = models.IntegerField(default = 1, choices = SNP_MODES)
    conMode = models.IntegerField(default = 1, choices = CON_MODE)
    scriptparams = models.TextField(blank = True, null = True)
    modifydate = models.DateTimeField(auto_now = True)
    single = models.BooleanField(default = False)
