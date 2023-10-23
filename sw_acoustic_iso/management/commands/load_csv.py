import csv
from itertools import islice
from django.conf import settings
from django.core.management.base import BaseCommand
from sw_acoustic_iso.models import Materials


class Command(BaseCommand):
    help = "Load data from csv file"
    
    def handle(self, *args, **kwargs):
        datafile = settings.BASE_DIR / 'data' / 'materiales_db.csv'
        
        with open(datafile, 'r') as csvfile:
            reader = csv.DictReader(islice(csvfile, 0, None))
            for row in reader:
                Materials.objects.get_or_create(material=row['Material'], 
                                                density=float(row['Densidad']),
                                                young_module=float(row['Modulo de Young']),
                                                loss_factor=float(row['Factor de perdidas']),
                                                poisson_module=float(row['Modulo Poisson']))