from cgbind import build
import numpy as np
from cgbind.linker import Linker
from cgbind.x_motifs import Xmotif


def test_get_fitted_linker_coords():

    l = Linker(arch_name='m2l4')

    # Generate a fictitious linker structure with two donor atoms..
    l.xyzs = [['N', 0.0, 0.0, 0.0], ['C', 5.0, 0.0, 0.0], ['N', 10.0, 0.0, 0.0]]

    l.x_motifs = [Xmotif(atom_ids=[0], coords=[np.array([0.0, 0.0, 0.0])]),
                  Xmotif(atom_ids=[2], coords=[np.array([10.0, 0.0, 0.0])])]

    # Template coords are aligned along the x axis
    template_coords = [np.array([0.0, 0.0, 0.0]), np.array([10.0, 0.0, 0.0])]

    coords_to_fit = [np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 10.0])]

    current_xyzs = [['X', 2.0, 0.0, 0.0]]

    ideal_linker_coords = np.array([[10., 0., 0.],
                                    [10., 0., 5.],
                                    [10., 0., 10.]])

    actual_transformed_coords = build.get_fitted_linker_coords(linker=l, template_x_coords=template_coords,
                                                               coords_to_fit=coords_to_fit, current_xyzs=current_xyzs)

    assert np.abs(np.sum(ideal_linker_coords - actual_transformed_coords)) < 0.001
