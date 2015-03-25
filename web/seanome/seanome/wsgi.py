"""
WSGI config for seanome project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.6/howto/deployment/wsgi/
"""

import os
import sys


project_path = os.path.split(os.path.split(os.path.realpath(__file__))[0])[0]
parent_path, project_name = os.path.split(project_path)

if(parent_path not in sys.path):
        sys.path.append(parent_path)

if(project_path not in sys.path):
        sys.path.append(project_path)

os.environ['DJANGO_SETTINGS_MODULE'] = '%s.settings'%(project_name)

from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()



#import os
#os.environ.setdefault("DJANGO_SETTINGS_MODULE", "seanome.settings")

#from django.core.wsgi import get_wsgi_application
#application = get_wsgi_application()


"""
WSGI config for seanome project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.6/howto/deployment/wsgi/
"""

