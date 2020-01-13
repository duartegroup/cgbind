from copy import deepcopy
import numpy as np
from scipy.spatial import distance_matrix
from cgbind.log import logger
from cgbind.geom import get_centered_matrix
from cgbind.geom import get_rot_mat_kabsch
from cgbind.geom import xyz2coord
from cgbind.x_motifs import get_shifted_template_x_motif_coords


def get_fitted_linker_coords_and_repulsion(linker, template_x_coords, coords_to_fit, current_xyzs):
    """
    For a linker get the best mapping onto a list of template X coords (e.g. NCN motifs in a pyridyl donor) these will
    this can be achieved in normal or reverse order of the coordinates as to maximise the distance to the rest of the
    metallocage structure. Also returns a measure of the repulsion to the rest of the cage structure

    :param linker: (Linker object)
    :param template_x_coords: (list(np.ndarray))
    :param coords_to_fit: (list(np.ndarray))
    :param current_xyzs: (list(list))
    :return: (list(np.ndarray)), (float)
    """
    curr_coords = xyz2coord(xyzs=current_xyzs)
    min_repulsion, best_linker_coords = 99999999999.9, None

    for coords in [coords_to_fit, list(reversed(coords_to_fit))]:
        new_linker_coords, _ = get_template_fitted_coords_and_cost(linker=linker,
                                                                   template_x_coords=template_x_coords,
                                                                   coords_to_fit=coords)
        if len(current_xyzs) == 0:
            return new_linker_coords, 0.0

        repulsion = np.sum(np.power(distance_matrix(new_linker_coords, curr_coords), -6))

        # Add the linker with the least repulsion to the rest of the structure
        if repulsion < min_repulsion:
            best_linker_coords = new_linker_coords
            min_repulsion = repulsion

    if best_linker_coords is None:
        logger.error('Fitted linker coords could not be found')

    return best_linker_coords, repulsion


def get_template_fitted_coords_and_cost(linker, template_x_coords, coords_to_fit, return_cost=False):
    """
    Get the coordinates of a linkers that are fitted to a template of X motifs

    :param linker: (Linker object)
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


def get_linker_xyzs_to_add_and_cost(linker, template_linker, curr_xyzs):
    """
    Get the xyzs of a linker that is fitted to a template_linker object and the associated cost function – i.e. the
    repulsion to the current cage structure

    :param linker: (Linker)
    :param template_linker: (Template.Linker)
    :param curr_xyzs: (list(list))
    :return: list(list)), float
    """

    # Ensure the shift amount dr is set
    if linker.dr is None:
        logger.error('Cannot build a cage dr was None')
        return None, 9999999.9

    # Copy the linker so no attributes are modified in place..
    new_linker = deepcopy(linker)

    # Expand the template by an about dr
    shifted_coords = get_shifted_template_x_motif_coords(linker_template=template_linker, dr=new_linker.dr)
    x_coords = [new_linker.coords[atom_id] for motif in new_linker.x_motifs for atom_id in motif.atom_ids]

    linker_coords, cost = get_fitted_linker_coords_and_repulsion(linker=new_linker,
                                                                 template_x_coords=shifted_coords,
                                                                 coords_to_fit=x_coords, current_xyzs=curr_xyzs)
    logger.info(f'Repulsive cost for adding the linker is {cost:.5f}')

    linker_xyzs = [[new_linker.xyzs[i][0]] + linker_coords[i].tolist() for i in range(new_linker.n_atoms)]
    return linker_xyzs, cost
