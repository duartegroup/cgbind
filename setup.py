from setuptools import setup
from Cython.Build import cythonize
from setuptools.extension import Extension
import sys
import os

here = os.path.abspath(os.path.dirname(__file__))

extensions = []

# By default don't build the ESP C++ extension on Windows
if '--build_ext' in sys.argv or sys.platform != 'win32':

    # Try to build the extension from the Cython generated C
    ext_folder = os.path.join(here, 'cgbind', 'ext')
    for file_name in ('esp_gen', 'conf_fit'):

        if os.path.exists(os.path.join(ext_folder, f'{file_name}.c')):
            ext = 'c'
        else:
            ext = 'pyx'

        extensions.append(Extension(file_name,
                                    [f'cgbind/ext/{file_name}.{ext}']))

if '--build_ext' in sys.argv:
    sys.argv.remove('--build_ext')


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
