from django.shortcuts import render
from sw_acoustic_iso.models import Materials
from sw_acoustic_iso.forms import MaterialsPanelForm, DimensionsPanel
from sw_acoustic_iso.acoustic import Panel
import plotly.graph_objects as go

# Create your views here.
def index(request):
    
    if request.method == "POST":
        material_form = MaterialsPanelForm(request.POST)
        dimensions_form = DimensionsPanel(request.POST)
        if material_form.is_valid() & dimensions_form.is_valid():
            
            material = Materials.objects.get(id=request.POST.get('material'))
            lx = request.POST.get('l_x')
            ly = request.POST.get('l_y')
            thickness = request.POST.get('thickness')
            
            panel = Panel(material.material, material.density, material.young_module, material.loss_factor, material.poisson_module, lx, ly, thickness)
            
            f_per_thirds = [20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000]
            
            f_cremer, r_cremer = panel.cremer_model(f_per_thirds)
            f_davy, r_davy = panel.davy_model() 
            f_sharp, r_sharp = panel.sharp_model(f_per_thirds)
            
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=f_cremer, y=r_cremer, name="Cremer"), )
            fig.add_trace(go.Scatter(x=f_davy, y=r_davy, name="Davy"))
            fig.add_trace(go.Scatter(x=f_sharp, y=r_sharp, name="Sharp"))
            # fig = px.line(x=f_cremer, y=r_cremer, log_x=True, color_discrete_sequence=['red'], labels={'x': 'Frecuencia [Hz]', 'y': 'R [dB]'})
            # fig.update_layout(legend_title_text='Aislamiento ', legend, autosize=True)
            # fig.update_xaxes(tickvals=f_per_thirds, ticktext=f_per_thirds)
            fig.update_xaxes(type="log", tickvals=f_per_thirds, ticktext=f_per_thirds)
            fig.update_layout(
                legend=dict(
                    title="Modelos",
                    orientation="h",
                    # entrywidth=70,
                    yanchor="bottom",
                    y=1.02,
                    xanchor="right",
                    x=1
                ),
                
                # legend=dict(
                #     x=0,
                #     y=1,
                #     title="Modelos",
                #     traceorder="reversed",
                # ),
                xaxis_title = "R [dB]",
                yaxis_title = "Frecuencia [Hz]",
                template = 'plotly_white'
            )
            
            fig_html = fig.to_html()

            
            return render(request, 'base.html', {'material_form': material_form, 'dimensions_form': dimensions_form, 'fig_html' : fig_html, 'stiffness' : panel.stiffness, 'mass_sup':panel.mass_sup, 'freq_res': panel.freq_res, 'freq_critic': panel.freq_critic, 'freq_density':panel.freq_density})
            
        else:
            print("Error en la validacion")
            #agregar mensaje de error
            pass
            # return render(request, 'base.html', {'material_form': material_form, 'dimensions_form': dimensions_form})
            
    elif request.method == "GET":
        dimensions_form = DimensionsPanel()
        material_form = MaterialsPanelForm()
        return render(request, 'base.html', {'material_form': material_form, 'dimensions_form': dimensions_form})
    else:
        return render(request, 'base.html')
        
    
    
