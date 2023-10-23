from django.shortcuts import render
from sw_acoustic_iso.models import Materials
from sw_acoustic_iso.forms import MaterialsPanelForm, DimensionsPanel

# Create your views here.
def index(request):
    
    material_form = MaterialsPanelForm()
    dimensions_form = DimensionsPanel()
    
    return render(request, 'base.html', {'material_form': material_form, 'dimensions_form': dimensions_form})