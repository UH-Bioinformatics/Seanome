# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding field 'job.buildmode'
        db.add_column(u'seanomeUI_job', 'buildmode',
                      self.gf('django.db.models.fields.IntegerField')(default=1),
                      keep_default=False)

        # Adding field 'job.repeatfilter'
        db.add_column(u'seanomeUI_job', 'repeatfilter',
                      self.gf('django.db.models.fields.IntegerField')(default=0),
                      keep_default=False)


    def backwards(self, orm):
        # Deleting field 'job.buildmode'
        db.delete_column(u'seanomeUI_job', 'buildmode')

        # Deleting field 'job.repeatfilter'
        db.delete_column(u'seanomeUI_job', 'repeatfilter')


    models = {
        u'seanomeUI.job': {
            'Meta': {'object_name': 'job'},
            'buildmode': ('django.db.models.fields.IntegerField', [], {'default': '1'}),
            'conMode': ('django.db.models.fields.IntegerField', [], {'default': '1'}),
            'data': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'interrupt': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'repeatfilter': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'searchMode': ('django.db.models.fields.IntegerField', [], {'default': '1'}),
            'sftp_id': ('django.db.models.fields.IntegerField', [], {'db_index': 'True'}),
            'snpcaller': ('django.db.models.fields.IntegerField', [], {'default': '1'}),
            'stage': ('django.db.models.fields.IntegerField', [], {'default': '0'})
        }
    }

    complete_apps = ['seanomeUI']