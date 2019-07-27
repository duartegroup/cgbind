import os
from cgbind import *
here = os.path.dirname(os.path.abspath(__file__))


def test_gen_cage():

    os.chdir(os.path.join(here, 'gen_cage'))
    cage_c1 = gen_cage(linker_name='L-1', linker_smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1')

    assert cage_c1.linker.charge == 0
    assert cage_c1.linker.name == 'L-1'
    assert len(cage_c1.linker.xyzs) == 34

    assert cage_c1.reasonable_geometry is True
    assert cage_c1.energy is None
    assert cage_c1.metal_charge == 2
    assert cage_c1.charge == 4
    assert cage_c1.metal == 'Pd'
    assert cage_c1.arch == M2L4
    assert cage_c1.name == 'cage_L-1'
    assert len(cage_c1.xyzs) == 138
