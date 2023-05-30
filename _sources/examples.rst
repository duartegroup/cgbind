Examples
========

.. currentmodule:: cgbind

A few examples for generating and analysing metallocages. The full examples are available on
`GitHub <https://github.com/duartegroup/cgbind/tree/master/examples>`_.

Generating a Metallocage
------------------------

To generate an archetype  Pd\ :sub:`2`\L\ :sub:`4` metallocage (`EZEVAI <https://dx.doi.org/10.5517/ccdc.csd.cc1m3h4c>`_).

.. code-block:: python

    >>> import cgbind
    >>> linker = cgbind.Linker(smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
    >>>                        arch_name='m2l4')
    >>> cage = cgbind.Cage(linker=linker, metal='Pd')

The :class:`Linker` holds the linker as a :class:`Molecule` and is required to build a :class:`Cage`. If a cage is
initialised from a single linker a homoleptic metallocage will be constructed. To generate a heteroleptic metallocage
the cage needs multple linkers, for example

    >>> alt_linker = cgbind.Linker(smiles='C1(C#CC2=NC(C#CC3=CC=CN=C3)=CC=C2)=CC=CN=C1',
    >>>                            arch_name='m2l4', name='alt_linker')
    >>> heteroleptic_cage = cgbind.Cage(linkers=[linker, linker, alt_linker, alt_linker],
    >>>                                 metal='Pd')

To generate a standard .xyz file that can be visualised with e.g. `Avogadro <https://avogadro.cc>`_ use

    >>> cage.print_xyz_file()



Obtaining Metallocage Properties
--------------------------------

To obtain properties of a specific :class:`Cage` e.g. the Pd-Pd distance, cavity volume (in Å\ :sup:`3`) and maximum
escape sphere (Å\ :sup:`3`) as a window size probe

    >>> m_m_dist = cage.get_m_m_dist()
    >>> vol = cage.get_cavity_vol()
    >>> escape_sphere = cage.get_max_escape_sphere()

To check all the current methods available

    >>> methods = [method for method in dir(cgbind.Cage) if callable(getattr(cgbind.Cage, method))]
    >>> print(methods)


M4L6 Metallocages
------------------

Two M\ :sub:`4`\L\ :sub:`6` templates are available in cgbind the first suitable for catechol derived linkers
(`arch_name='m4l6'`), for example the cage published by `Raymond <http://dx.doi.org/10.1021/ja411631v>`_ . Also
available is a template suitable for bipyridyl based linkers (`arch_name='m4l6n'`, also referred to as M4L6'),
for example published by `Nitchke <https://doi.org/10.1021/ja402084x>`_. To generate a bipyridyl-based
Fe(II)\ :sub:`4`\L\ :sub:`6` metallocage over two cores:

.. code-block:: python

    >>> from cgbind import Linker, Cage, Config
    >>> Config.n_cores = 2
    >>> linker = Linker(smiles='C1(C2=CC=C(C#CC3=CN=C(C4=NC=CC=C4)C=C3)C=N2)=NC=CC=C1', arch_name='m4l6n')
    >>> cage = Cage(linker=linker, metal='Fe', metal_charge=2)


Due to the conformational flexibility of some linkers a larger number of conformers may be required to find one
with the correct geometry. To request more than the default 300 conformers initialise a linker with
`Linker(..., n_confs=300)`. Also possible is the ETKDGv2 algorithm in RDKit does not generate a conformer that affords a
sensible structure, to use an earlier version `Linker(..., use_etdg_confs=True)`.


Adding Architectures
--------------------

Extending *cgbind* to other architectures is possible using the add_template.py
`script <https://github.com/duartegroup/cgbind/blob/master/scripts/add_template.py>`_.
For example, the M24L48 metallocage from `Fujita <https://doi.org/10.1126/science.1188605>`_ can
be added as a template by downloading the crystal structure deposited in the
`CCDC <https://www.ccdc.cam.ac.uk/structures/Search?Ccdcid=765717&DatabaseToSearch=Published>`_
using Open → CSD entry in external viewer to download the structure as a .mol2. Then::

    $ python /path/to/cgbind/scripts/add_template.py WUQROQ.mol2 --name m24l48

Creating a M24L48 cage is then as easy as:

.. code-block:: python

    >>> from cgbind import Linker, Cage
    >>> linker = Linker(smiles='n1ccc(cc1)-c2sc(cc2)-c3ccncc3', arch_name='m24l48')
    >>> cage = Cage(linker=linker, metal='Pd')
    >>> cage.print_xyz_file()





