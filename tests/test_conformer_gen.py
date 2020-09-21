from cgbind import Linker, Cage
from cgbind.linker_conformers import set_best_fit_linker
from cgbind.input_output import xyz_file_to_atoms
from conf_fit import rotate_coords
import networkx as nx
import numpy as np
import os


def test_molecule_split():

    linker = Linker(
                    # smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                    # smiles='C1=CN=CC(C#CC2=CC=CC(C#CC3=CN=CC=C3)=C2)=C1',
                    # smiles='C1(C#CC2=CC=CN=C2)=CC(C#CC3=CN=CC=C3)=CC=C1',
                    # smiles='CC(C)C1=CC=C(/C=N/C2=CC=C(C3=CC=C(C(C=CC(C4=CC=C(/N=C/C5=CC=C(C(C)C)C=N5)C=C4)=C6)=C6C7)C7=C3)C=C2)N=C1',
                    # smiles='CC([C@H]1C2)(C)[C@@H](C1)C3=C2C=CC(C(C=C4)=NC=C4C5=CC=C(C6=CC=C(C[C@@H]7C(C)(C)[C@H]8C7)C8=N6)N=C5)=N3',
                    # smiles='O=C(NC1=CC=CC2=C1C=CC=C2NC(C3=C([O-])C([O-])=CC=C3)=O)C4=C([O-])C([O-])=CC=C4',
                    smiles='C1(CCCCCCCC2=CC=CN=C2)=CN=CC=C1',
                    arch_name='m2l4',
                    n_confs=1)
    linker.print_xyz_file(filename='unrot.xyz')
    linker.set_ranked_linker_possibilities(metal='Pd')
    linker.possibilities[0].print_xyz_file(filename='rot.xyz')

    cage = Cage(linker, metal='Fe', metal_charge=2)
    cage.print_xyz_file()


def test_cython_rotation():

    xyz_path = os.path.join('tests', 'data', 'ch_linker.xyz')
    linker = Linker(arch_name='m2l4',
                    filename=xyz_path)

    coords = linker.get_coords()
    linker.graph.remove_edge(2, 3)
    components = list(nx.connected_components(linker.graph))
    linker.graph.add_edge(2, 3)

    rot_coords = rotate_coords(coords,
                               (coords[2] - coords[3]),
                               np.pi,
                               coords[3],
                               list(components[0]))

    linker.set_atoms(coords=rot_coords)
    linker.print_xyz_file(filename='tmp.xyz')
