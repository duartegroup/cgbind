from cgbind.molecule import Molecule
import numpy as np
import os
from cgbind import xtb

here = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()


def test_molecule():

    methane = Molecule(name='methane', smiles='C')
    assert methane.charge == 0
    assert methane.mult == 1
    assert methane.solvent is None
    assert methane.n_atoms == 5
    assert methane.reasonable_geometry is True
    assert methane.energy is None

    coords = methane.get_coords()
    assert coords.shape == (5, 3)

    # Centering the molecule should make the centroid at the origin
    methane.centre()
    coords = methane.get_coords()
    assert np.linalg.norm(np.average(coords, axis=0)) < 1E-1

    # Reset the molecules atoms from a altered set of coordinates
    coords[0] = np.ones(3)
    methane.set_atoms(coords=coords)
    assert all(methane.atoms[0].coord == np.ones(3))

    methane.print_xyz_file()
    assert os.path.exists('methane.xyz')


def test_molecule_sp_opt():
    os.chdir(os.path.join(here, 'data'))

    methane = Molecule(name='methane', smiles='C')

    methane.singlepoint(method=xtb)
    assert methane.energy == -4.173842879099

    methane.optimise(method=xtb)
    assert methane.energy == -4.175218496669

    os.chdir(cwd)


def test_rdkit_props():
    os.chdir(os.path.join(here, 'data'))

    methane = Molecule(name='methane', smiles='C')
    assert methane.n_rot_bonds == 0
    assert methane.n_h_donors == 0
    assert methane.n_h_acceptors == 0
    assert len(methane.bonds) == 4
    assert (0, 1) in methane.bonds
    assert methane.mol_obj is not None

    charges = methane.get_charges(estimate=True)
    assert len(charges) == methane.n_atoms

    # Should be no large charges in methane
    assert all(-1.0 < c < 1.0 for c in charges)

    # Run an XTB calculation to get the charges
    charges = methane.get_charges()
    assert len(charges) == methane.n_atoms
    os.chdir(cwd)


def test_confs():

    butane = Molecule(name='butane', smiles='CCCC', n_confs=5)
    assert len(butane.conformers) > 1
    assert butane.conformers[0].n_atoms == butane.n_atoms


def test_molecule_from_file():

    methane = Molecule(name='methane', filename=os.path.join(here, 'data', 'methane.xyz'))
    assert methane.n_atoms == 5
    assert methane.charge == 0
    assert methane.mult == 1

    methane = Molecule(name='methane', filename=os.path.join(here, 'data', 'methane.mol2'))
    assert methane.n_atoms == 5

    methane = Molecule(name='methane', filename=os.path.join(here, 'data', 'methane.mol'))
    assert methane.n_atoms == 5