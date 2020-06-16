from setuptools import setup
from Cython.Build import cythonize
from setuptools.extension import Extension
import sys
import os

here = os.path.abspath(os.path.dirname(__file__))

extensions = []

# By default don't build the ESP C++ extension on Windows
if '--esp_gen' in sys.argv or sys.platform != 'win32':

    # Try to build the extension from the Cython generated C
    if os.path.exists(os.path.join(here, 'cgbind', 'ext', 'esp_gen.c')):
        ext = 'c'
    else:
        ext = 'pyx'

    extensions.append(Extension('esp_gen', [f'cgbind/ext/esp_gen.{ext}']))

if '--esp_gen' in sys.argv:
    sys.argv.remove('--esp_gen')


setup(name='cgbind',
      version='1.0.0a',
      description='Metallocage construction and binding affinity calculations',
      packages=['cgbind'],
      package_data={'': ['lib/*']},
      ext_modules=cythonize(extensions, language_level="3", annotate=True),
      url='https://github.com/duartegroup/cgbind',
      license='MIT',
      author='Tom Young',
      author_email='tom.young@chem.ox.ac.uk')
