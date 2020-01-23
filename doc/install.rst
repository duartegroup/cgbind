Install
=======

cgbind is built in Python 3 and has been tested with v. 3.7. An `anaconda <https://www.anaconda.com/distribution>`_ or
`miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ installation is highly recommended to satisfy the
dependencies easily.

.. note::
    The `autode <https://duartegroup.github.io/autodE/install.html>`_ package is a dependency and must be installed
    prior to ``cgbind``, but does not require any electronic structure theory packages.


To install ``cgbind`` generate with a conda environment::

    $ conda config --append channels conda-forge
    $ conda create -n autode_env --file requirements.txt
    $ conda activate cgbind_env

then, in the cloned repository (i.e. ``git clone https://github.com/duartegroup/cgbind.git`` and ``cd cgbind``)::

    $ python setup.py install

