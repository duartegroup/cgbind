from setuptools import setup

setup(
    name='cgbind',
    version='1.0.0',
    packages=['cgbind'],
    include_package_data=True,
    package_data={'': ['lib/*.obj']},
    url='',
    license='MIT',
    author='Tom Young',
    author_email='tom.young@chem.ox.ac.uk',
    description='Cage Binding Affinity Calculations'
)
