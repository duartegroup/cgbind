from copy import deepcopy
import numpy as np
from cgbind.geom import get_centered_matrix
from cgbind.geom import get_rot_mat_kabsch
from cgbind.geom import xyz2coord


def get_template_fitted_coords_and_cost(linker, template_x_coords, coords_to_fit):
    """
    Get the coordinates of a linkers that are fitted to a template of XEEX motifs
    :param linker: (object)
    :param template_x_coords: (list(np.ndarray))
    :param coords_to_fit: (list(np.ndarray)) must have len() = len(linker_template.x_xyzs)
    :return: (np.ndarray) n_atoms x 3
    """

    # Construct the P matrix in the Kabsch algorithm
    p_mat = deepcopy(coords_to_fit)
    p_centroid = np.average(p_mat, axis=0)
    p_mat_trans = get_centered_matrix(p_mat)

    # Construct the P matrix in the Kabsch algorithm
    q_mat = deepcopy(template_x_coords)
    q_centroid = np.average(q_mat, axis=0)
    q_mat_trans = get_centered_matrix(q_mat)

    # Get the optimum rotation matrix
    rot_mat = get_rot_mat_kabsch(p_mat_trans, q_mat_trans)

    # Apply to get the new set of coordinates
    new_linker_coords = np.array([np.matmul(rot_mat, coord - p_centroid) + q_centroid
                                  for coord in xyz2coord(linker.xyzs)])

    # Compute the cost function = (r - r_ideal)^2
    x_atom_ids = [x for x_motif in linker.x_motifs for x in x_motif]
    new_p_mat = np.array([new_linker_coords[i] for i in x_atom_ids])
    cost = np.sum(np.square(new_p_mat - q_mat))

    return new_linker_coords, cost


def get_cost_xyzs(r, cage, linker, template, coords_to_fit, return_xyzs=False):
    """
    For an  metallocage compute the cost function, which is the sum of the square differences between both the
    template XEEX motifs and the fitted linkers. This cost will be minimised with respect to the size of the
    which is quantified by the centroid-M distance

    :param r: (float)
    :param cage: (object)
    :param linker: (object)
    :param template: (object)
    :param coords_to_fit: (list(np.ndarray)) list of coords to fit onto the template
    :param return_xyzs:
    :return:
    """

    # Modify a copy of the template inplace with an r which corresponds to the centroid-M distance
    mod_template = deepcopy(template)
    mod_template.set_geom(r=r)

    # Add the xyzs of the metals to the list of xyzs, inserting the appropriate metal atom label
    xyzs = [[cage.metal] + m_template_xyzs[1:] for m_template_xyzs in mod_template.Metals.xyzs]

    cost = 0.0
    for n, ligand_template in enumerate(mod_template.Ligands):
        ligand_coords, l_cost = get_template_fitted_coords_and_cost(linker=linker, linker_template=ligand_template,
                                                                    coords_to_fit=coords_to_fit)
        xyzs += [[linker.xyzs[i][0]] + ligand_coords[i].tolist() for i in range(linker.n_atoms)]
        cost += l_cost

        if n == 0:
            from cgbind.input_output import xyzs2xyzfile
            xyzs2xyzfile(xyzs=xyzs, basename='tmp4')

    if return_xyzs:
        return cost, xyzs
    else:
        return cost
