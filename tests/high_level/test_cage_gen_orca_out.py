import os
from cgbind import *
here = os.path.dirname(os.path.abspath(__file__))


def test_gen_cage_orca_out():

    os.chdir(os.path.join(here, 'gen_cage'))
    Config.code = 'orca'
    Config.path_to_orca = '/Users/tom/opt/orca_4_1/orca'

    cage_c1 = gen_cage(linker_name='L-1', linker_smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                       opt_linker=False, opt_cage=False, metal_label='Pd',
                       metal_charge=2, n_cores_pp=1, sp_cage=True)

    assert cage_c1.energy == -3770.064764747303
    assert len(cage_c1.xyzs) == 138
