Examples
========

.. currentmodule:: cgbind

A few examples for generating and analysing metallocages. The full examples are available on
`GitHub <https://github.com/duartegroup/cgbind/tree/master/examples>`_.

Generating a Metallocage
------------------------

To generate a simple  Pd\ :sub:`2` L\ :sub:`4` metallocage

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

To generate a standard .xyz file that can be visualised with `Avogadro <https://avogadro.cc>`_ use

    >>> cage.print_xyzfile()


Obtaining Metallocage Properties
--------------------------------

To obtain properties of a specific :class:`Cage` e.g. the Pd-Pd distance, cavity volume (in Å\ :sup:`3`) and maximum
escape sphere (Å\ :sup:`3`) as a window size probe

    >>> m_m_dist = cage.get_m_m_dist()
    >>> vol = cage.get_cavity_vol()
    >>> escape_sphere = cage.get_max_escape_sphere()

To check all the current methods available

    >>> methods = [method for method in dir(cage) if callable(getattr(cgbind.Cage, method))]
    >>> print(methods)
