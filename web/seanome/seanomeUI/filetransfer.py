import os
import urllib
#from django.http import HttpResponse
from django.utils.encoding import smart_str

from zipstream import *
from django.http import StreamingHttpResponse


def zip_path(path, archive_dir_name = ""):
    """Recursively generate data to add directory tree or file pointed to by                 
    path to the archive. Results in archive containing                                       
    archive_dir_name/basename(path)                                                          
    archive_dir_name/basename(path)/*                                                        
    archive_dir_name/basename(path)/*/*                                                      
    .                                                                                        
    .                                                                                        
    .                                                                                        
    path -- path to file or directory                                                        
    archive_dir_name -- name of containing directory in archive                              
    """
    if os.path.isdir(path):
        dir_name = os.path.basename(path)
        for name in os.listdir(path):
            r_path = os.path.join(path, name)
            r_archive_dir_name = os.path.join(archive_dir_name, dir_name)
            for f,p in zip_path(r_path, r_archive_dir_name):
                yield (f,p)
    else:
        archive_path = os.path.join(archive_dir_name, os.path.basename(path))
        yield (path, archive_path)



def streamZipArchive(path, name):

    zstream = ZipFile(fileobj = None, compression = ZIP_DEFLATED)
    for f, p in zip_path(path):
        zstream.write(f,p)
    response = StreamingHttpResponse(zstream, mimetype='application/zip')
    response['Content-Disposition'] = 'attachment; filename="%s"' %(smart_str(name) )
    response['Expires'] = 0
    response['Accept-Ranges'] = 'bytes'
    response['Cache-Control'] = "private"
    response['Pragma'] = 'private'
    return response


def streamArchive(zstream, name):
    response = StreamingHttpResponse(zstream, mimetype='application/zip')
    response['Content-Disposition'] = 'attachment; filename="%s"' %(smart_str(name) )
    response['Expires'] = 0
    response['Accept-Ranges'] = 'bytes'
    response['Cache-Control'] = "private"
    response['Pragma'] = 'private'
    return response
