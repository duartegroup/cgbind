from setuptools import setup, Command
from Cython.Build import cythonize
from setuptools.extension import Extension
import sys

# By default don't build the ESP C++ extension on Windows
if '--esp_gen' not in sys.argv and sys.platform == 'win32':
    extensions = []
else:
    extensions = [Extension('esp_gen', ['cgbind/ext/esp_gen.pyx'])]

if '--esp_gen' in sys.argv:
    sys.argv.remove('--esp_gen')


setup(name='cgbind',
      version='1.0.0',
      packages=['cgbind'],
      include_package_data=True,
      package_data={'': ['lib/*']},
      ext_modules=cythonize(extensions, language_level="3", annotate=True),
      url='',
      license='MIT',
      author='Tom Young',
      author_email='tom.young@chem.ox.ac.uk')

