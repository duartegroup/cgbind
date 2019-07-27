from cgbind import *


def test_gen_linker_no_smiles():

    linker = Linker(name='linker')
    assert linker.smiles is None
    assert linker.name == 'linker'
    assert linker.charge == 0
    assert linker.xyzs is None
    assert linker.energy is None
    assert linker.conf_ids is None
