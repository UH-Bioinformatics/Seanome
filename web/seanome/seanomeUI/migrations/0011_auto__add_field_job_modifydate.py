# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding field 'job.modifydate'
        db.add_column(u'seanomeUI_job', 'modifydate',
                      self.gf('django.db.models.fields.DateTimeField')(auto_now=True, default=datetime.datetime(2014, 11, 6, 0, 0), blank=True),
                      keep_default=False)


    def backwards(self, orm):
        # Deleting field 'job.modifydate'
        db.delete_column(u'seanomeUI_job', 'modifydate')


    models = {
        u'seanomeUI.job': {
            'Meta': {'object_name': 'job'},
            'conMode': ('django.db.models.fields.IntegerField', [], {'default': '1'}),
            'data': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'interrupt': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'modifydate': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'scriptparams': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'searchMode': ('django.db.models.fields.IntegerField', [], {'default': '1'}),
            'sftp_id': ('django.db.models.fields.IntegerField', [], {'db_index': 'True'}),
            'snpcaller': ('django.db.models.fields.IntegerField', [], {'default': '1'}),
            'stage': ('django.db.models.fields.IntegerField', [], {'default': '0'})
        }
    }

    complete_apps = ['seanomeUI']