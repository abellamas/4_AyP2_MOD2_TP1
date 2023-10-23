from django.db import models
from django.contrib.auth.models import AbstractUser

# Create your models here.
class Materials(models.Model):
    material = models.CharField(verbose_name="Material", max_length=100)
    density = models.IntegerField(verbose_name="Densidad")
    young_module = models.FloatField(verbose_name="Módulo de Young")
    loss_factor = models.FloatField(verbose_name="Factor de perdidas")
    poisson_module = models.FloatField(verbose_name="Módulo de Poisson")
    
    def __str__(self):
        return self.material
    
    
class MaterialsPanel(models.Model):
    material = models.ForeignKey(Materials, on_delete=models.SET_NULL, blank=True, null=True)