# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding field 'job.searchMode'
        db.add_column(u'seanomeUI_job', 'searchMode',
                      self.gf('django.db.models.fields.IntegerField')(default=1),
                      keep_default=False)


    def backwards(self, orm):
        # Deleting field 'job.searchMode'
        db.delete_column(u'seanomeUI_job', 'searchMode')


    models = {
        u'seanomeUI.job': {
            'Meta': {'object_name': 'job'},
            'data': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'interrupt': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'searchMode': ('django.db.models.fields.IntegerField', [], {'default': '1'}),
            'sftp_id': ('django.db.models.fields.IntegerField', [], {'db_index': 'True'}),
            'stage': ('django.db.models.fields.IntegerField', [], {'default': '0'})
        }
    }

    complete_apps = ['seanomeUI']