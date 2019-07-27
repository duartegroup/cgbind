import os
from cgbind import *
here = os.path.dirname(os.path.abspath(__file__))


def test_gen_cage_bad_linker():

    os.chdir(os.path.join(here, 'gen_cage'))
    cage_x = gen_cage(linker_name='LX', linker_smiles='C1=CC=CN=C1')

    assert cage_x.reasonable_geometry is False
    assert cage_x.xyzs is None
