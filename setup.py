from setuptools import setup
from Cython.Build import cythonize
from setuptools.extension import Extension
import sys
import os

here = os.path.abspath(os.path.dirname(__file__))

extensions = []

# By default don't build the ESP C++ extension on Windows
if '--esp_gen' in sys.argv or sys.platform != 'win32':

    # Generate with cython always (to keep up with cython versions)
    extensions.append(Extension('esp_gen', [f'cgbind/ext/esp_gen.pyx']))

if '--esp_gen' in sys.argv:
    sys.argv.remove('--esp_gen')


setup(name='cgbind',
      version='1.0.3',
      description='Metallocage construction and binding affinity calculations',
      packages=['cgbind'],
      package_data={'': ['lib/*']},
      ext_modules=cythonize(extensions, language_level="3", annotate=True),
      url='https://github.com/duartegroup/cgbind',
      license='MIT',
      author='Tom Young',
      author_email='tom.young@chem.ox.ac.uk',
      install_requires=['Cython'],
      python_requires=">3.6")
