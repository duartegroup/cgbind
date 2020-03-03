Install
=======

cgbind is built in Python 3 and has been tested with v. 3.7. An `anaconda <https://www.anaconda.com/distribution>`_ or
`miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ installation is highly recommended to easily install the
dependencies.

.. note::
    The `autode <https://duartegroup.github.io/autodE/install.html>`_ package is a dependency and must be installed
    prior to ``cgbind``, but does not require any electronic structure theory packages to also be installed.

To install ``cgbind`` inside a conda environment first close the repository and ``cd`` there::

    $ git clone https://github.com/duartegroup/cgbind.git
    $ cd cgbind

then create and activate a virtual environment with the appropriate dependencies::

    $ conda config --append channels conda-forge
    $ conda create -n cgbind_env --file requirements.txt
    $ conda activate cgbind_env

finally::

    $ python setup.py install

