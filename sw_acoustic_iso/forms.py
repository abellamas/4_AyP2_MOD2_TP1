from django import forms
from django.contrib.auth.forms import UserCreationForm
from sw_acoustic_iso.models import Materials, MaterialsPanel


class MaterialsPanelForm(forms.ModelForm):
    class Meta:
        model = MaterialsPanel
        fields = ('material',)
        
class DimensionsPanel(forms.Form):
    l_x = forms.FloatField(required=True, min_value=0, label="Largo x")
    l_y = forms.FloatField(required=True, min_value=0, label="Largo y")
    thickness = forms.FloatField(required=True, min_value=0, label="Espesor")