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
    assert cage_c1.cage_substrate_name is None
    assert cage_c1.cage_substrate_xyzs is None
    assert cage_c1.cage_substrate_energy is None
    assert cage_c1.charge == 4
    assert cage_c1.metal == 'Pd'
    assert cage_c1.arch == 'm2l4'
    assert cage_c1.name == 'cage_L-1'
    assert len(cage_c1.xyzs) == 138
    assert Config.n_cores == 1


def test_gen_cage_bad_linker():

    os.chdir(os.path.join(here, 'gen_cage'))

    cage_x = gen_cage(linker_name='LX', linker_smiles='C1=CC=CN=C1')

    assert cage_x.linker.xyzs is None
    assert cage_x.reasonable_geometry is False
    assert cage_x.xyzs is None


# def test_gen_m4l6_cage():
#
#     os.chdir(os.path.join(here, 'gen_cage'))
#
#     cage_m4l6 = gen_cage(linker_name='L_m4l6',
#                          linker_smiles='O=C(C1=C([O-])C([O-])=CC=C1)NC2=CC=CC3=C2C=CC=C3NC(C4=C([O-])C([O-])=CC=C4)=O',
#                          metal_label='Ga',
#                          metal_charge=3,
#                          arch='m4l6')
#
#     assert cage_m4l6.linker.charge == -4
#
#     assert cage_m4l6.charge == -12
#     assert cage_m4l6.metal == 'Ga'
#     assert cage_m4l6.arch == 'm4l6'
#     assert cage_m4l6.energy is None
#     assert len(cage_m4l6.xyzs) == 280


def test_gen_cage_orca_out():

    os.chdir(os.path.join(here, 'gen_cage'))
    Config.code = 'orca'
    Config.path_to_orca = '/Users/tom/opt/orca_4_1/orca'

    cage_c1 = gen_cage(linker_name='L-1', linker_smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                       opt_linker=False, opt_cage=False, metal_label='Pd',
                       metal_charge=2, n_cores_pp=1, sp_cage=True)

    assert cage_c1.energy == -3770.064764747303
    assert len(cage_c1.xyzs) == 138


def test_n_cores_pp():

    from cgbind.parallel import calc_n_cores_pp

    dict1, dict2 = {'test1': 'none', 'test2': 'none'}, {'test1': 'none', 'test2': 'none'}

    # Default total number of cores is 1
    assert calc_n_cores_pp(dict1, dict2) == 1

    Config.n_cores = 2
    assert calc_n_cores_pp(dict1) == 1
    assert calc_n_cores_pp(dict1, dict2) == 1

    Config.n_cores = 4
    assert calc_n_cores_pp(dict1) == 2
    assert calc_n_cores_pp(dict1, dict2) == 1


def test_gen_m2l4_cage_subst_complex_bq():

    os.chdir(os.path.join(here, 'gen_cage_subst'))

    cage_c1_bq, bq = gen_cage_subst_complex(linker_name='L-1',
                                            linker_smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                                            substrate_name='bq', substrate_smiles='O=C1C=CC(C=C1)=O')

    assert cage_c1_bq.reasonable_geometry is True
    assert cage_c1_bq.energy is None
    assert len(cage_c1_bq.cage_substrate_xyzs) == 150
    assert cage_c1_bq.cage_substrate_energy is None
    assert cage_c1_bq.charge == 4
    assert cage_c1_bq.metal == 'Pd'
    assert cage_c1_bq.name == 'cage_L-1'
    assert cage_c1_bq.cage_substrate_name == 'cage_L-1_bq'
    assert cage_c1_bq.cage_subst_reasonable_geometry is True
    assert len(cage_c1_bq.xyzs) == 138

    assert bq.name == 'bq'
    assert len(bq.xyzs) == 12
    assert bq.n_heteroatoms == 2
    assert bq.x_x_atom_ids[0] is not None
    assert bq.x_x_atom_ids[1] is not None
    assert bq.energy is None


def test_gen_m2l4_cage_subst_complex_acetone():

    os.chdir(os.path.join(here, 'gen_cage_subst'))

    cage_c1_bq, acetone = gen_cage_subst_complex(linker_name='L-1',
                                                 linker_smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                                                 substrate_name='acetone', substrate_smiles='CC(C)=O')

    assert cage_c1_bq.reasonable_geometry is True
    assert cage_c1_bq.energy is None
    assert len(cage_c1_bq.cage_substrate_xyzs) == 138 + 10
    assert cage_c1_bq.cage_substrate_name == 'cage_L-1_acetone'
    assert cage_c1_bq.cage_subst_reasonable_geometry is True

    assert acetone.name == 'acetone'
    assert len(acetone.xyzs) == 10
    assert acetone.n_heteroatoms == 1
    assert acetone.x_atom_ids[0] is not None
    assert acetone.energy is None


def test_gen_m2l4_cage_subst_complex_methane():

    os.chdir(os.path.join(here, 'gen_cage_subst'))

    cage_c1_bq, methane = gen_cage_subst_complex(linker_name='L-1',
                                                 linker_smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                                                 substrate_name='methane', substrate_smiles='C')

    assert cage_c1_bq.reasonable_geometry is True
    assert cage_c1_bq.energy is None
    assert len(cage_c1_bq.cage_substrate_xyzs) == 138 + 5
    assert cage_c1_bq.cage_substrate_name == 'cage_L-1_methane'
    assert cage_c1_bq.cage_subst_reasonable_geometry is True

    assert len(methane.xyzs) == 5
    assert methane.n_heteroatoms == 0
    assert methane.energy is None


def test_gen_linker_no_smiles():

    linker = Linker(name='linker')
    assert linker.smiles is None
    assert linker.name == 'linker'
    assert linker.charge == 0
    assert linker.xyzs is None
    assert linker.energy is None
    assert hasattr(linker, 'conf_ids') is False


def test_gen_substrate_no_smiles():

    substrate = Substrate(name='substrate')
    assert substrate.name == 'substrate'
    assert substrate.charge == 0
    assert substrate.energy is None
    assert substrate.xyzs is None
    assert substrate.x_x_atom_ids is None
    assert substrate.x_atom_ids is None

