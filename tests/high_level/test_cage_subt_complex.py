import os
from cgbind import *
here = os.path.dirname(os.path.abspath(__file__))


def test_gen_m2l4_cage_subst_complex_bq():

    os.chdir(os.path.join(here, 'gen_cage_subst'))

    cage_c1_bq = gen_cage_subst_complex(linker_name='L-1', linker_smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                                        substrate_name='bq', substrate_smiles='O=C1C=CC(C=C1)=O')

    assert cage_c1_bq.cage.reasonable_geometry is True
    assert cage_c1_bq.cage.energy is None
    assert len(cage_c1_bq.xyzs) == 150
    assert cage_c1_bq.energy is None
    assert cage_c1_bq.charge == 4
    assert cage_c1_bq.cage.metal == 'Pd'
    assert cage_c1_bq.cage.name == 'cage_L-1'
    assert cage_c1_bq.name == 'cage_L-1_bq'
    assert cage_c1_bq.reasonable_geometry is True
    assert len(cage_c1_bq.cage.xyzs) == 138

    assert cage_c1_bq.substrate.name == 'bq'
    assert len(cage_c1_bq.substrate.xyzs) == 12
    assert cage_c1_bq.substrate.n_heteroatoms == 2
    assert cage_c1_bq.substrate.x_x_atom_ids[0] is not None
    assert cage_c1_bq.substrate.x_x_atom_ids[1] is not None
    assert cage_c1_bq.substrate.energy is None


def test_gen_m2l4_cage_subst_complex_acetone():

    os.chdir(os.path.join(here, 'gen_cage_subst'))

    cage_c1_acetone = gen_cage_subst_complex(linker_name='L-1',
                                             linker_smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                                             substrate_name='acetone', substrate_smiles='CC(C)=O')

    assert cage_c1_acetone.cage.reasonable_geometry is True
    assert cage_c1_acetone.energy is None
    assert len(cage_c1_acetone.xyzs) == 138 + 10
    assert cage_c1_acetone.name == 'cage_L-1_acetone'
    assert cage_c1_acetone.reasonable_geometry is True

    assert cage_c1_acetone.substrate.name == 'acetone'
    assert len(cage_c1_acetone.substrate.xyzs) == 10
    assert cage_c1_acetone.substrate.n_heteroatoms == 1
    assert cage_c1_acetone.substrate.x_atom_ids[0] is not None
    assert cage_c1_acetone.substrate.energy is None


def test_gen_m2l4_cage_subst_complex_methane():

    os.chdir(os.path.join(here, 'gen_cage_subst'))

    cage_c1_methane = gen_cage_subst_complex(linker_name='L-1',
                                             linker_smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                                             substrate_name='methane', substrate_smiles='C')

    assert cage_c1_methane.cage.reasonable_geometry is True
    assert cage_c1_methane.cage.energy is None
    assert len(cage_c1_methane.xyzs) == 138 + 5
    assert cage_c1_methane.name == 'cage_L-1_methane'
    assert cage_c1_methane.reasonable_geometry is True

    assert len(cage_c1_methane.substrate.xyzs) == 5
    assert cage_c1_methane.substrate.n_heteroatoms == 0
    assert cage_c1_methane.substrate.energy is None
