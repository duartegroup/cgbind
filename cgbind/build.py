from copy import deepcopy
import numpy as np
from cgbind.geom import get_centered_matrix
from cgbind.geom import get_rot_mat_kabsch
from cgbind.geom import xyz2coord


def get_template_fitted_coords_and_cost(linker, template_x_coords, coords_to_fit, return_cost=False):
    """
    Get the coordinates of a linkers that are fitted to a template of X motifs

    :param linker: (object)
    :param template_x_coords: (list(np.ndarray))
    :param coords_to_fit: (list(np.ndarray)) must have len() = len(linker_template.x_xyzs)
    :param return_cost: (bool) return just the cost function, which is the sum of squares of ∆dists
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

    if return_cost:
        new_p_mat = np.array([np.matmul(rot_mat, coord) for coord in p_mat_trans])
        cost = np.sum(np.square(np.array([np.linalg.norm(new_p_mat[i] - q_mat_trans[i]) for i in range(len(coords_to_fit))])))
        return cost

    # Apply to get the new set of coordinates
    new_linker_coords = np.array([np.matmul(rot_mat, coord - p_centroid) + q_centroid
                                  for coord in xyz2coord(linker.xyzs)])

    # Compute the cost function = (r - r_ideal)^2
    x_atom_ids = [x for x_motif in linker.x_motifs for x in x_motif.atom_ids]
    new_p_mat = np.array([new_linker_coords[i] for i in x_atom_ids])
    cost = np.sum(np.square(np.array([np.linalg.norm(new_p_mat[i] - q_mat[i]) for i in range(len(coords_to_fit))])))

    return new_linker_coords, cost
