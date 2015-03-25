from django import forms
from .models import *
from django_select2.widgets import Select2Widget
from simplemathcaptcha.fields import MathCaptchaField
import json



class JobStarter(forms.Form):
    def __init__(self, *args, **kwargs):
        super(JobStarter, self).__init__(*args, **kwargs)
        self.fields['captcha'] = MathCaptchaField()


class minimumCluster(forms.Form):
    def __init__(self, maxsize, *args, **kwargs):
        super(minimumCluster, self).__init__(*args, **kwargs)
        self.fields['minimum'] = forms.IntegerField(label = "Minimum number of instances", min_value = 0, max_value = maxsize, default = 0)

