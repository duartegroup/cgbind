from cgbind import atoms


def test_atoms():

    assert atoms.get_atomic_number(atom_label='P') == 15
    # For unknown atoms types the default atomic number is 6
    assert atoms.get_atomic_number(atom_label='XX') == 6

    assert 11.99 < atoms.get_atomic_mass(atom_label='C') < 12.02
    # For unknown atoms types the default atomic mass is 10
    assert 9 < atoms.get_atomic_mass(atom_label='XX') < 11

    assert 1.1 < atoms.get_vdw_radii(atom_label='H') < 1.3
    # For unknown atoms types the default van der Walls radii is 1.5 Å
    assert 1.4 < atoms.get_vdw_radii(atom_label='XX') < 1.6
