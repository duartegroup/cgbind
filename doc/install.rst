Install
=======

cgbind is built in Python 3 and has been tested with v. 3.7. An `anaconda <https://www.anaconda.com/distribution>`_ or
`miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ installation is highly recommended to easily install the
dependencies.

.. note::
    The `autode <https://duartegroup.github.io/autodE/install.html>`_ package is an optional dependency and must be
    installed prior to ``cgbind``, to perform semi-emprical, DFT and ab initio calculations.

Mac OSX / Linux
---------------

To install ``cgbind`` inside a conda environment first clone the repository and ``cd`` there::

    $ git clone https://github.com/duartegroup/cgbind.git
    $ cd cgbind


then install the appropriate dependencies (you may want to create a new virtual environment)::

    $ conda config --append channels conda-forge
    $ conda install --file requirements.txt

finally::

    $ python setup.py install



Windows
--------

On Windows without a ``git`` installation ``cgbind`` can be installed with `anaconda <https://www.anaconda.com/distribution>`_
by, on the GitHub `page <https://github.com/duartegroup/cgbind>`_ using Clone or download → Download ZIP then
extracting it. Then, open an anaconda command prompt and ``cd`` to the directory and proceed as above e.g.::

    $ cd Downloads/cgbind-master/
    $ conda config --append channels conda-forge
    $ conda install --file requirements.txt
    $ python setup.py install
