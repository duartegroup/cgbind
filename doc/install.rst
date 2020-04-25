Install
=======

cgbind is built in Python 3 and has been tested with v. 3.7. An `anaconda <https://www.anaconda.com/distribution>`_ or
`miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ installation is highly recommended to easily install the
dependencies.

.. note::
    The `autode <https://duartegroup.github.io/autodE/install.html>`_ package is an optional dependency to perform
    semi-emprical, DFT and ab initio calculations.

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

The above commands assume you have extracted the zip to ``C:\Users\yourusername\Downloads``.

.. note::
    By default Windows doesn't have a C++ compiler required to build the extension to generate electrostatic potentials,
    as such this feature is not enabled by default. If a C++ compiler is available install with ``python setup.py install --esp_gen`` to build the ESP extension.

