from django.conf import settings
from django.conf.urls import patterns, include, url
from django.contrib.auth.views import login, logout, password_change, password_change_done


urlpatterns = patterns('seanomeUI.views',
                       #url(r'^$', 'staticPage', {'template':'index.html'}, name ="home"),
                       url(r'^$', 'beginJob', {'template':'jobstarter.html'}, name ="startjob"),
                       url(r'^run/job/(?P<jid>\d+)/$', 'runJob', {'template':'jobRunner.html'}, name ="runjob"),
                       url(r'^results/job/(?P<jid>\d+)/$', 'completedJob', {'template':'jobComplete.html'}, name ="completedjob"),
                       url(r'^results/job/filter/summary/(?P<jid>\d+)/(?P<count>\d+)/$', 'resultFilter', {'template':'filterResults.html'}, name ="filterResults"),
                       url(r'^results/job/filter/alignment/(?P<jid>\d+)/(?P<count>\d+)/$', 'viewAlignments', {'template':'viewalignment.html'}, name ="viewaln"),
                       url(r'^results/job/filter/alignment/(?P<jid>\d+)/(?P<count>\d+)/2/$', 'viewAlignments', {'template':'viewalignment2.html'}, name ="viewaln2"),

                       url(r'^download/(?P<jid>\d+)/(?P<count>\d+)/(?P<tid>\d+)/$', 'MSADownload', name ="msa_download"),
                       url(r'^download/SNP/(?P<jid>\d+)/(?P<count>\d+)/$', 'downloadvcf', name ="snp_download"),
                       url(r'^cancel/job/(?P<jid>\d+)/$', 'cancelJob', {'template':'jobRunner.html'}, name ="canceljob"),
                       url(r'^ajax/refresh/directory/(?P<jid>\d+)/$', 'ajaxRefreshFileDir', name ="ajax_dir_refresh"),
                       url(r'^ajax/job/status/(?P<jid>\d+)/$', 'ajaxJobStatus', name ="ajax_status"),
                       url(r'^ajax/msa/(?P<jid>\d+)/(?P<ident>\d+)/(?P<clean>[01])/$', 'ajaxGetMSA', name ="ajax_msa"),
                       url(r'^ajax/snp/(?P<jid>\d+)/(?P<ident>\d+)/$', 'ajaxGetSNPs', name ="ajax_snps"),
)

#urlpatterns += patterns('',
#)

if settings.DEBUG:
    import debug_toolbar
    urlpatterns += patterns('',
                        url(r'^__debug__/', include(debug_toolbar.urls)),
                        )
