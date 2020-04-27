from cgbind import cage
from cgbind.linker import Linker
import numpy as np
import os

here = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()


def test_init_cage():

    linker = Linker(arch_name='m2l4', filename=os.path.join(here, 'data', 'ch_linker.xyz'), name='tmp')

    c = cage.Cage(linker=linker, metal='Pd', metal_charge=2)
    assert c._is_linker_reasonable(linker=linker)

    assert np.abs(np.sum(c.get_centroid() - np.array([0.25825, 17.5488,  12.52415]))) < 0.001

    c.print_xyz_file()
    assert os.path.exists('cage_tmp.xyz')
    os.remove('cage_tmp.xyz')

    assert c.get_metal_atom_ids() == [136, 137]
    assert 110 < c.get_cavity_vol() < 135       # Å^3
    assert 11.5 < c.get_m_m_dist() < 12.5       # Å

    # The below are none as the structure is not initialised with SMILES string, thus there are no RDKit objects
    assert c.get_num_rot_bonds() is None
    assert c.get_num_h_bond_donors() is None
    assert c.get_num_h_bond_acceptors() is None

    assert 0.0 < c.get_max_escape_sphere() < 1000.0
