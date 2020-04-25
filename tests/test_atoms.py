from cgbind import atoms
import numpy as np


def test_atom_class():

    x = atoms.Atom('X', 0.0, 0.0, 0.0)
    assert x.coord.shape == (3,)
    assert all(x.coord == np.zeros(3))
    assert len(str(x).split()) == 4     # Atom symbol, x, y, z


def test_atoms():

    assert atoms.get_atomic_number(atom=atoms.Atom('P', 0.0, 0.0, 0.0)) == 15
    # For unknown atoms types the default atomic number is 6
    assert atoms.get_atomic_number(atom=atoms.Atom('XX', 0.0, 0.0, 0.0)) == 6

    assert 11.99 < atoms.get_atomic_mass(atom=atoms.Atom('C', 0.0, 0.0, 0.0)) < 12.02
    # For unknown atoms types the default atomic mass is 10
    assert 9 < atoms.get_atomic_mass(atom=atoms.Atom('XX', 0.0, 0.0, 0.0)) < 11

    assert 1.1 < atoms.get_vdw_radii(atom=atoms.Atom('H', 0.0, 0.0, 0.0)) < 1.3
    # For unknown atoms types the default van der Walls radii is 1.5 Å
    assert 1.4 < atoms.get_vdw_radii(atom=atoms.Atom('XX', 0.0, 0.0, 0.0)) < 1.6

    assert atoms.get_max_valency(atom=atoms.Atom('H', 0.0, 0.0, 0.0)) == 1
    assert atoms.get_max_valency(atom=atoms.Atom('C', 0.0, 0.0, 0.0)) == 4
    assert atoms.get_max_valency(atom=atoms.Atom('XX', 0.0, 0.0, 0.0)) == 6
