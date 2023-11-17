from django import forms
from django.contrib.auth.forms import UserCreationForm
from sw_acoustic_iso.models import Materials, MaterialsPanel


class MaterialsPanelForm(forms.ModelForm):
    class Meta:
        model = MaterialsPanel
        fields = ('material',)
        widgets = {
            'material': forms.Select(attrs={'class':"form-select"}),
        }
        
class DimensionsPanel(forms.Form):
    l_x = forms.FloatField(required=True, min_value=0, label="Largo x [m] ")
    l_y = forms.FloatField(required=True, min_value=0, label="Largo y [m] ")
    thickness = forms.FloatField(required=True, min_value=0, label="Espesor [m]")
    
    def __init__(self, *args, **kwargs):
        super(DimensionsPanel, self).__init__(*args, **kwargs)
        for visible in self.visible_fields():
            visible.field.widget = forms.TextInput()
            visible.field.widget.attrs['class'] = 'form-control'