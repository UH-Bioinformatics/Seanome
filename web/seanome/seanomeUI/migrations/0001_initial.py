# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'job'
        db.create_table(u'seanomeUI_job', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('sftp_id', self.gf('django.db.models.fields.IntegerField')(db_index=True)),
        ))
        db.send_create_signal(u'seanomeUI', ['job'])


    def backwards(self, orm):
        # Deleting model 'job'
        db.delete_table(u'seanomeUI_job')


    models = {
        u'seanomeUI.job': {
            'Meta': {'object_name': 'job'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'sftp_id': ('django.db.models.fields.IntegerField', [], {'db_index': 'True'})
        }
    }

    complete_apps = ['seanomeUI']