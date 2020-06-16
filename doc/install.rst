Install
=======

cgbind is built in Python 3 and has been tested with v. 3.7. An `anaconda <https://www.anaconda.com/distribution>`_ or
`miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ installation is highly recommended to easily install the
dependencies.

.. note::
    The `autode <https://duartegroup.github.io/autodE/install.html>`_ package is an optional dependency to perform
    semi-emprical, DFT and ab initio calculations.


PIP
---

If the requirements (rdkit, numpy, scipy, networkx, Cython) are already available ``cgbind`` can be install with
`pip <https://pypi.org/project/pip/>`_::

    $ pip install cgbind



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


The whole process including installing miniconda is shown `here <https://youtu.be/R-J6vJeydAE>`_.

Windows
--------

On Windows without a ``git`` installation ``cgbind`` can be installed with `anaconda <https://www.anaconda.com/distribution>`_
by, on the GitHub `page <https://github.com/duartegroup/cgbind>`_ using Clone or download → Download ZIP then
extracting it. Then, open an anaconda command prompt and ``cd`` to the directory and proceed as above e.g.::

    $ cd Downloads/cgbind-master/
    $ conda config --append channels conda-forge
    $ conda install --file requirements.txt
    $ python setup.py install

The above commands assume you have extracted the zip to ``C:\Users\yourusername\Downloads``.

.. note::
    By default Windows doesn't have a C++ compiler required to build the extension to generate electrostatic potentials,
    as such this feature is not enabled by default. If a C++ compiler is available install with ``python setup.py install --esp_gen`` to build the ESP extension.


Common Problems
---------------

Conda is not found::

    $ conda config --append channels conda-forge
    conda command not found

First, make sure either miniconda or anaconda installed then ensure you've closed and reopened your terminal after
installation.

A required module is not found (setuptools, rdkit, Cython)::

    $ python setup.py install
    Traceback (most recent call last):
      File "setup.py", line 2, in <module>
        from Cython.Build import cythonize
    ImportError: No module named Cython.Build

Check that you've run ``conda install --file requirements.txt`` and that you're using the conda version of python (i.e.
running ``which python`` returns /some/path/(ana/mini)conda/bin/python). If ``which python`` returns /usr/bin/python or
another non-conda version of python then you may need to activate the base conda environment with ``conda activate base``.
