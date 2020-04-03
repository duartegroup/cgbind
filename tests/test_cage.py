from cgbind import cage
from cgbind.linker import Linker
import numpy as np
import os


def test_init_cage():

    xyzs = [['C', -7.9509, -0.8834, -1.2393],
            ['C', -6.5711, -0.5321, -0.8275],
            ['N', -5.5405, -1.3889, -1.0461],
            ['C', -4.2735, -1.0577, -0.6639],
            ['C', -6.3397, 0.6785, -0.2171],
            ['C', -5.0552, 1.005, 0.1659],
            ['C', -3.9914, 0.1453, -0.048],
            ['C', -2.6533, 0.4763, 0.3467],
            ['C', -1.5326, 0.7749, 0.6762],
            ['C', -0.1859, 1.1592, 1.0685],
            ['C', -0.0459, 2.1746, 2.0116],
            ['C', 1.2143, 2.57, 2.4118],
            ['C', 2.2826, 1.9175, 1.8375],
            ['C', 2.1355, 0.9134, 0.9047],
            ['C', 3.2504, 0.2331, 0.3059],
            ['C', 4.2005, -0.3287, -0.1852],
            ['C', 5.342, -0.9885, -0.7497],
            ['C', 6.618, -0.4662, -0.5961],
            ['C', 7.7239, -1.095, -1.1364],
            ['C', 7.5747, -2.2669, -1.8461],
            ['N', 6.3448, -2.766, -1.9913],
            ['C', 5.2637, -2.1619, -1.4698],
            ['N', 0.8931, 0.5524, 0.5367],
            ['H', -8.6072, -1.0446, -0.3666],
            ['H', -7.9507, -1.732, -1.9501],
            ['H', -8.3884, -0.0185, -1.8025],
            ['H', -3.4907, -1.7798, -0.8635],
            ['H', -7.1718, 1.3484, -0.051],
            ['H', -4.8819, 1.9542, 0.6431],
            ['H', -0.952, 2.6333, 2.4111],
            ['H', 1.3044, 3.3543, 3.1398],
            ['H', 3.2855, 2.1893, 2.1163],
            ['H', 6.7024, 0.4617, -0.0297],
            ['H', 8.7211, -0.6823, -1.0136],
            ['H', 8.461, -2.7635, -2.2745],
            ['H', 4.2649, -2.5855, -1.5992]]

    l = Linker(arch_name='m2l4', xyzs=xyzs, name='tmp')

    c = cage.Cage(linker=l, metal='Pd', metal_charge=2)
    assert c._is_linker_reasonable(linker=l)

    assert np.abs(np.sum(c.get_centroid() - np.array([ 0.25825, 17.5488,  12.52415]))) < 0.001

    c.print_xyzfile()
    assert os.path.exists('cage_tmp.xyz')
    os.remove('cage_tmp.xyz')

    assert c.get_metal_atom_ids() == [144, 145]
    assert 122 < c.get_cavity_vol() < 135       # Å^3
    assert 11.5 < c.get_m_m_dist() < 12.5       # Å

    # The below are none as the structure is not initialised with SMILES string, thus there are no RDKit objects
    assert c.get_num_rot_bonds() is None
    assert c.get_num_h_bond_donors() is None
    assert c.get_num_h_bond_acceptors() is None

    assert 0.0 < c.get_max_escape_sphere() < 1000.0
