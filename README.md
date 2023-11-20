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
5. Carga de materiales de .xlsx -> .csv a base de datos db.sqlite3

Los materiales a cargar desde un excel .xlsx deben transformarse a un .csv delimitado por comas que deberá estar en el siguiente directorio */data/materiales_db.csv*. 
Debe tenerse cuidado con los carácteres utilizados, por ej. evitar el uso de ñ, diéresis, acentos, etc.  Tambien se recomienda revisar que todos los números transformados de excel .xlsx a .csv coincidan, en especial las unidades pequeñas.
Para ejecutar la carga del .csv a la base de datos por defecto db.sqlite3 se debe ejecutar el comando:

```console
python manage.py load_csv
```
Una vez cargado no debe entregar ningun mensaje por consola, y deben realizarse las migraciones para aplicar los cambios del ORM a la base de datos.

```console
python manage.py makemigrations
```
```console
python manage.py migrate
```
6. Exportación de archivo Excel con resultados

Para exportar los resultados se usa un template base que se encuentra en */data/template.xlsx*.
Por el momento, en cada simulación se genera internamente un archivo excel */data/templates.xlsx* que no se descarga hasta que el usuario presiona el botón de _*Exportar*_.
Dicho archivo se descarga en el directorio predeterminado por el navegador a utilizar y se guardará con el siguiente formato: *aislamiento_ruido_HH MM SS.xlsx*.

Alerta: dado que el archivo excel a descargar siempre se genera, si se presiona exportar antes de realizar la simulación puede que se obtenga la última simulación realizada.

# Posibles mejoras a futuro

1. Permitir la carga y modificación de materiales desde la vista en navegador.
2. Evitar la generación del archivo excel siempre que se simule y solo hacerlo cuando se presione Exportar, incluyendo las caracteristicas del material en el nombre del archivo.
3. Correcciones de la simulación física en el codigo que se encuentren con las revisiones y uso.