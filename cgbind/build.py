from copy import deepcopy
import numpy as np
from scipy.optimize import minimize, Bounds
from scipy.spatial import distance_matrix
from cgbind.log import logger
from cgbind.geom import get_centered_matrix
from cgbind.geom import get_rot_mat_kabsch
from cgbind.geom import xyz2coord
from cgbind.x_motifs import get_shifted_template_x_motif_coords


def get_fitted_linker_coords_and_cost(linker, template_x_coords, coords_to_fit, current_xyzs):
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
    min_cost, best_linker_coords = 99999999999.9, None

    for coords in [coords_to_fit, list(reversed(coords_to_fit))]:
        new_linker_coords, cost = get_kfitted_coords_and_cost(linker=linker,
                                                              template_x_coords=template_x_coords,
                                                              coords_to_fit=coords)
        if len(current_xyzs) == 0:
            return new_linker_coords, 0.0

        repulsion = np.sum(np.power(distance_matrix(new_linker_coords, curr_coords), -12))

        # Add the linker with the least repulsion to the rest of the structure
        if repulsion + cost < min_cost:
            best_linker_coords = new_linker_coords
            min_cost = repulsion + cost

    if best_linker_coords is None:
        logger.error('Fitted linker coords could not be found')

    return best_linker_coords, min_cost


def get_kfitted_coords_and_cost(linker, template_x_coords, coords_to_fit, return_cost=False):
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

    linker_coords, cost = get_fitted_linker_coords_and_cost(linker=new_linker, template_x_coords=shifted_coords,
                                                            coords_to_fit=x_coords, current_xyzs=curr_xyzs)
    logger.info(f'Repulsive + fitting cost for adding the linker is {cost:.5f}')

    if linker_coords is None:
        logger.error('Linker coords were None')
        return None, cost

    linker_xyzs = [[new_linker.xyzs[i][0]] + linker_coords[i].tolist() for i in range(new_linker.n_atoms)]
    return linker_xyzs, cost


def cost_fitted_x_motifs(dr, linker, linker_template, x_coords):
    """
    For a linker compute the cost function (RMSD) for fitting all the coordinates in the x motifs to a template which
    which be shifted by dr in the corresponding shift_vec

    :param linker: (object)
    :param linker_template: (object)
    :param x_coords: (list(np.ndarray))
    :param dr: (float)
    :return:
    """

    shifted_coords = get_shifted_template_x_motif_coords(linker_template=linker_template, dr=dr)
    cost = get_kfitted_coords_and_cost(linker, template_x_coords=shifted_coords, coords_to_fit=x_coords,
                                       return_cost=True)
    return cost


def get_new_linker_and_cost(linker_xyzs, linker, x_motifs, template_linker):
    """
    For a set of linker xyzs e.g. one conformer and a linker object together with a list of x motifs in the
    structure return cost function associated with fitting the X motifs to the template

    :param linker_xyzs: (list(list))
    :param linker: (Linker)
    :param x_motifs: (list(Xmotif))
    :param template_linker: (Template.Linker)
    :return: (tuple(Linker, float)) New linker and cost
    """

    coords = xyz2coord(linker_xyzs)
    x_coords = [coords[atom_id] for motif in x_motifs for atom_id in motif.atom_ids]

    # Minimise the cost function as a function of dr in Å
    min_result = minimize(cost_fitted_x_motifs, x0=np.array([1.0]), args=(linker, template_linker, x_coords),
                          method='L-BFGS-B', tol=1e-3, bounds=Bounds(-100.0, 10.0))

    # Set attributes of the new linker
    new_linker = deepcopy(linker)
    new_linker.dr = min_result.x[0]
    new_linker.xyzs = linker_xyzs
    new_linker.coords = coords
    new_linker.centroid = np.average(new_linker.coords, axis=0)
    new_linker.x_motifs = x_motifs

    return new_linker, min_result.fun


def build_homoleptic_cage(cage, max_cost):
    """
    Construct the geometry (xyzs) of a homoleptic cage

    :param cage: (Cage)
    :param max_cost: (float) Maximum cost to break out of the loop over
    :return:
    """

    # Get the list of Linkers ordered by the best fit to the template
    linker_conf_list = cage.linkers[0].get_ranked_linker_conformers(metal=cage.metal)
    logger.info(f'Have {len(linker_conf_list)} conformers to fit')

    min_cost, best_linker = 99999999.9, None
    xyzs = []

    # TODO: Parallelise
    for linker in linker_conf_list:
        cage_cost = 0.0

        for i, template_linker in enumerate(cage.cage_template.linkers):

            linker_xyzs, cost = get_linker_xyzs_to_add_and_cost(linker, template_linker, curr_xyzs=xyzs)
            cage_cost += cost

            if linker_xyzs is None:
                logger.error('Failed to get linker')
                break

            xyzs += linker_xyzs

        if cage_cost < min_cost:
            min_cost = cage_cost
            best_linker = deepcopy(linker)

        if cage_cost < max_cost:
            logger.info(f'Total L-L repulsion + fit to template in building cage is {cage_cost:.2f}')
            break

        else:
            xyzs = []

    if len(xyzs) == 0:
        if best_linker is None:
            logger.error('Could not achieve the required cost threshold for building the cage')
            return None
        else:
            logger.warning('Failed to reach the threshold. Returning the cage that minimises the L-L repulsion')

    cage.dr = best_linker.dr
    xyzs = []
    for i, template_linker in enumerate(cage.cage_template.linkers):
        linker_xyzs, _ = get_linker_xyzs_to_add_and_cost(best_linker, template_linker, curr_xyzs=xyzs)
        xyzs += linker_xyzs

    # Add the metals from the template shifted by dr
    for metal in cage.cage_template.metals:
        metal_coord = cage.dr * metal.shift_vec / np.linalg.norm(metal.shift_vec) + metal.coord
        xyzs.append([cage.metal] + metal_coord.tolist())

    cage.set_xyzs(xyzs)
    return None
