# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding field 'job.stage'
        db.add_column(u'seanomeUI_job', 'stage',
                      self.gf('django.db.models.fields.IntegerField')(default=0),
                      keep_default=False)


    def backwards(self, orm):
        # Deleting field 'job.stage'
        db.delete_column(u'seanomeUI_job', 'stage')


    models = {
        u'seanomeUI.job': {
            'Meta': {'object_name': 'job'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'sftp_id': ('django.db.models.fields.IntegerField', [], {'db_index': 'True'}),
            'stage': ('django.db.models.fields.IntegerField', [], {'default': '0'})
        }
    }

    complete_apps = ['seanomeUI']