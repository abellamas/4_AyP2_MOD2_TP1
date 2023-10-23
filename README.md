# Instalación

0. Abrir la terminal en el directorio donde se descarguen el repositorio.

1. Crear un entorno virtual, yo lo llamo '.venv' pero pueden ponerle cualquier nombre.
   
```console
python -m virtualenv .venv
```
2. Activar el entorno virtual
```console
.venv/scripts/activate
```
3. Verificar estar en la consola ubicados a la altura de la carpeta principal (donde está el archivo `manage.py`) e instalar todas las dependencias en el entorno virtual
```console
pip install -r requirements.txt
```
4. Iniciar el servidor
```console
python manage.py runserver
```


