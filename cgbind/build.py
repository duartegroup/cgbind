from copy import deepcopy
import numpy as np
from scipy.optimize import minimize, Bounds
from scipy.spatial import distance_matrix
from cgbind.log import logger
from cgbind.atoms import Atom
from cgbind.config import Config
from cgbind.geom import get_centered_matrix
from cgbind.geom import get_rot_mat_kabsch
from cgbind.x_motifs import get_shifted_template_x_motif_coords
from cgbind.molecule import BaseStruct


def get_fitted_linker_coords_and_cost(linker, template_x_coords, coords_to_fit, curr_coords):
    """
    For a linker get the best mapping onto a list of template X coords (e.g. NCN motifs in a pyridyl donor) these will
    this can be achieved in normal or reverse order of the coordinates as to maximise the distance to the rest of the
    metallocage structure. Also returns a measure of the repulsion to the rest of the cage structure

    :param linker: (Linker object)
    :param template_x_coords: (list(np.ndarray))
    :param coords_to_fit: (list(np.ndarray))
    :param curr_coords: (list(list))
    :return: (list(np.ndarray)), (float)
    """
    min_cost, best_linker_coords = 99999999999.9, None

    coord_set = [coords_to_fit]

    if Config.allow_permutations:
        logger.info('Allowing permutations of coords - only reverse')
        coord_set.append(list(reversed(coords_to_fit)))

    for coords in coord_set:
        new_linker_coords, cost = get_kfitted_coords_and_cost(linker=linker,
                                                              template_x_coords=template_x_coords,
                                                              coords_to_fit=coords)
        if len(curr_coords) == 0:
            return new_linker_coords, 0.0

        repulsion = np.sum(np.power(distance_matrix(new_linker_coords, curr_coords), -12))

        # Add the linker with the least repulsion to the rest of the structure
        if repulsion + cost < min_cost:
            best_linker_coords = new_linker_coords
            min_cost = repulsion + cost

    if best_linker_coords is None:
        logger.warning('Fitted linker coords could not be found')
        best_linker_coords = linker.get_coords()

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
    assert len(template_x_coords) == len(coords_to_fit)

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
    new_linker_coords = np.array([np.matmul(rot_mat, coord - p_centroid) + q_centroid for coord in linker.get_coords()])

    # Compute the cost function = (r - r_ideal)^2
    x_atom_ids = [x for x_motif in linker.x_motifs for x in x_motif.atom_ids]
    new_p_mat = np.array([new_linker_coords[i] for i in x_atom_ids])
    cost = np.sum(np.square(np.array([np.linalg.norm(new_p_mat[i] - q_mat[i]) for i in range(len(coords_to_fit))])))

    return new_linker_coords, cost


def get_linker_atoms_and_cost(linker, template_linker, current_atoms, x_coords=None):
    """
    Get the xyzs of a linker that is fitted to a template_linker object and the associated cost function – i.e. the
    repulsion to the current cage structure

    :param linker: (Linker)
    :param template_linker: (Template.Linker)
    :param curr_coords: (list(list))
    :return: list(list)), float
    """

    if x_coords is None:
        x_coords = linker.get_xmotif_coordinates()

    # Ensure the shift amount dr is set
    if linker.dr is None:
        logger.error('Cannot build a cage dr was None')
        return linker.atoms, 9999999.9

    # Expand the template by an about dr
    shifted_coords = get_shifted_template_x_motif_coords(linker_template=template_linker,
                                                         dr=linker.dr)

    linker_coords, cost = get_fitted_linker_coords_and_cost(linker=linker,
                                                            template_x_coords=shifted_coords,
                                                            coords_to_fit=x_coords,
                                                            curr_coords=[atom.coord for atom in current_atoms])
    logger.info(f'Repulsive + fitting cost for adding the linker is {cost:.5f}')

    atoms = [Atom(linker.atoms[i].label, coord=linker_coords[i])
             for i in range(linker.n_atoms)]

    return atoms, cost


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


def build_homoleptic_cage(cage, max_cost):
    """
    Construct the geometry (atoms) of a homoleptic cage

    :param cage: (Cage)
    :param max_cost: (float) Maximum cost to break out of the loop over
    :return:
    """

    # Get the list of Linkers ordered by the best fit to the template
    cage.linkers[0].set_ranked_linker_possibilities(metal=cage.metal)
    logger.info(f'Have {len(cage.linkers[0].possibilities)} linkers to fit')

    min_cost, best_linker = 99999999.9, None

    # For all the possible linker conformer / Xmotif set possibilities
    for linker in cage.linkers[0].possibilities:

        # Atoms for and cost in building this cage
        atoms, cage_cost = [], 0.0

        # Coordinates of the X motif atoms in this linker - used to rotate
        x_coords = linker.get_xmotif_coordinates()

        for i, template_linker in enumerate(cage.cage_template.linkers):
            linker_atoms, cost = get_linker_atoms_and_cost(linker,
                                                           template_linker,
                                                           atoms,
                                                           x_coords)
            cage_cost += cost
            atoms += linker_atoms

        if cage_cost < min_cost:
            min_cost = cage_cost
            best_linker = deepcopy(linker)

        if cage_cost < max_cost:
            logger.info(f'Total L-L repulsion + fit to template in building '
                        f'cage is {cage_cost:.2f}')
            break

    # If there is no break due to a small repulsion then build the best
    # possible cage
    if cage_cost > max_cost:
        if best_linker is None:
            logger.error('Could not achieve the required cost threshold for '
                         'building the cage')
            return None
        else:
            logger.warning('Failed to reach the threshold. Returning the cage '
                           'that minimises the L-L repulsion')
            atoms = []
            for i, template_linker in enumerate(cage.cage_template.linkers):
                linker_atoms, _ = get_linker_atoms_and_cost(best_linker,
                                                            template_linker,
                                                            atoms)
                atoms += linker_atoms

    # Set the delta r for the whole cage
    cage.dr = best_linker.dr

    # Add the metals from the template shifted by dr
    for metal in cage.cage_template.metals:
        metal_coord = cage.dr * metal.shift_vec / np.linalg.norm(metal.shift_vec) + metal.coord
        atoms.append(Atom(cage.metal, coord=metal_coord))

    cage.set_atoms(atoms)
    return None


def build_heteroleptic_cage(cage, max_cost):
    logger.info('Building a heteroleptic cage')
    logger.warning('Due to the very large space that needs to be minimised '
                   'only the *best* linker conformer is used')

    added_linkers, atoms = [], []

    for i, linker in enumerate(cage.linkers):

        linker.set_ranked_linker_possibilities(metal=cage.metal)

        linker_atoms, cost = get_linker_atoms_and_cost(linker.possibilities[0],
                                                       cage.cage_template.linkers[i],
                                                       atoms)

        logger.info(f'L-L repulsion + fit to template in building cage is {cost:.2f}')
        atoms += linker_atoms
        linker.dr = linker.possibilities[0].dr

    logger.warning('Heteroleptic cages will have the average dr of all linkers '
                   '- using the average')
    cage.dr = np.average(np.array([linker.dr for linker in cage.linkers]))

    # Add the metals from the template shifted by dr
    for metal in cage.cage_template.metals:
        metal_coord = cage.dr * metal.shift_vec / np.linalg.norm(metal.shift_vec) + metal.coord
        atoms.append(Atom(cage.metal, x=metal_coord[0], y=metal_coord[1], z=metal_coord[2]))

    cage.set_atoms(atoms)
    return None
