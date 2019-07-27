from cgbind import *


def test_gen_substrate_no_smiles():

    substrate = Substrate(name='substrate')
    assert substrate.name == 'substrate'
    assert substrate.charge == 0
    assert substrate.energy is None
    assert substrate.xyzs is None
    assert substrate.x_x_atom_ids is None
    assert substrate.x_atom_ids is None
