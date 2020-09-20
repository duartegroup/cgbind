from cgbind import Linker
from cgbind.linker_conformers import set_best_fit_linker
from cgbind.input_output import xyz_file_to_atoms
from conf_fit import rotate_coords
import networkx as nx
import numpy as np
import os


def test_molecule_split():

    linker = Linker(
                    # smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                    smiles='CCCC1=CN=CC(C#CC2=CC=CC(C#CC3=CN=CC=C3)=C2)=C1',
                    arch_name='m2l4',
                    use_fragment_conf=False,
                    n_confs=1)
    # linker.set_atoms(atoms=xyz_file_to_atoms('tests/data/ch_linker.xyz'))

    set_best_fit_linker(linker, linker.x_motifs,
                        template_linker=linker.cage_template.linkers[0])


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
