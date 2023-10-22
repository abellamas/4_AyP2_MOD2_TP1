from django.shortcuts import render
from sw_acoustic_iso.models import Materials

# Create your views here.
def index(request):
    
    materials = Materials.objects.all()
    
    return render(request, 'base.html')