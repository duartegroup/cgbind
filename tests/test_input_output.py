from cgbind import input_output
import os


def test_xyz_file():

    path_to_data = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))
    os.chdir(path_to_data)

    atoms = input_output.xyzfile_to_atoms(filename='cage_tmp.xyz')
    assert len(atoms) == 146


def test_mol_file():

    atoms = input_output.molfile_to_xyzs(filename='methane.mol')
    assert len(atoms) == 5


def test_mol2_file():

    atoms = input_output.mol2file_to_atoms(filename='methane.mol2')
    assert len(atoms) == 5
