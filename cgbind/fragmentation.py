from cgbind.geom import get_rot_mat_kabsch
from cgbind.geom import rotation_matrix
from cgbind.bonds import get_avg_bond_length
from cgbind.log import logger
from cgbind.utils import fast_xtb_opt
from scipy.spatial import distance_matrix
from scipy.optimize import minimize
from itertools import product
from copy import deepcopy
from time import time
import networkx as nx
import numpy as np


def get_frag_minimised_linker(linker, fragments, template_linker, x_motifs):
    """
    From a linker with a set of fragments, fit the fragments containing
    x-motifs using a defined dr to the template, then minimise the remaining
    fragments under a

    V = bonds + repulsion

    potential with respect to rotation and translation of the fragment.

    :param linker: (cgbind.linker.Linker)
    :param fragments: (list(list(int))) List of atom indexes of each fragment
    :param template_linker: (cgbind.templates.Linker)
    :param x_motifs: (list(cgbind.x_motifs.Xmotif))
    :return: (float or np.ndarray or None)
    """
    start_time = time()
    coords = np.array(linker.get_coords(), copy=True)

    fitted_idxs = []
    shift_vecs = []

    for i, fragment in enumerate(fragments):

        # If the fragment contains an x_motif this function fits it to the
        # template and returns the shift vector
        shift_vec = x_motif_vec_and_fit(linker, fragment, template_linker,
                                        coords, x_motifs)
        if shift_vec is not None:
            fitted_idxs.append(i)
            shift_vecs.append(shift_vec)

    # Fitted coordinates that are fixed up to the shift vector value
    x_frag_idxs = [np.array(fragment, dtype=np.int)
                   for i, fragment in enumerate(fragments) if i in fitted_idxs]

    flat_x_frag_idxs = np.array([idx for i in fitted_idxs for idx in fragments[i]],
                                dtype=np.int)

    # Deal with the remaining, not fitted fragment
    fragment = [frag for i, frag in enumerate(fragments)
                if i not in fitted_idxs][0]

    # Minimise the energy with respect to the remaining fragments
    fragment_idxs = np.array(fragment, dtype=np.int)
    ij_bonds = get_cross_bonds(linker, fragment, flat_x_frag_idxs)

    min_opt = None

    for _ in range(15):
        opt = minimize(rotated_fragment_linker_energy,
                       x0=np.concatenate((np.random.uniform(size=8, low=-1.0, high=2.0)),
                                         axis=None),
                       args=(coords,
                             fragment_idxs,
                             x_frag_idxs,
                             flat_x_frag_idxs,
                             ij_bonds,
                             shift_vecs),
                       method='BFGS')
        if min_opt is None or opt.fun < min_opt.fun:
            min_opt = opt

    # print(opt)
    energy = rotated_fragment_linker_energy(min_opt.x,
                                            coords,
                                            fragment_idxs,
                                            x_frag_idxs,
                                            flat_x_frag_idxs,
                                            ij_bonds,
                                            shift_vecs,
                                            modify=True)
    linker = deepcopy(linker)
    linker.dr = min_opt.x[0]
    linker.x_motifs = x_motifs
    linker.cost = energy
    linker.set_atoms(coords=coords)

    # Try to optimise slightly with XTB
    fast_xtb_opt(linker, n_cycles=10)

    logger.info(f'Generated linker in {time() - start_time:.3f} s')
    return linker


def rotated_fragment_linker_energy(x, coords, idxs_to_shift,
                                   x_motifs_idxs,
                                   flat_x_motifs_idxs,
                                   ij_bonds,
                                   shift_vecs,
                                   modify=False):
    """
    Rotate a fragment of the linker using a numpy array defining the shift
    vector (length 3) and the rotation (length 4)

    :param x:
    :param coords:
    :param idxs_to_shift:
    :param x_motifs_idxs:
    :param flat_x_motifs_idxs: 1D version of x_motifs_idxs
    :param shift_vecs: (list(np.ndarray)) Shift vectors for the first dimension
                       in x_motifs_idxs
    :param ij_bonds:
    :param modify: (bool) Modify the coordinates in place
    :return:
    """
    if not modify:
        coords = np.array(coords, copy=True)

    # Shift all the fragments with x motifs in them with dr with the the
    # associated shift vector

    dr = x[0]
    for i, x_motif_idxs in enumerate(x_motifs_idxs):
        coords[x_motif_idxs] += dr * shift_vecs[i]

    # Translate the shiftable portion
    coords[idxs_to_shift] += x[1:4]

    # And rotate it
    rot_mat = rotation_matrix(axis=np.concatenate((x[4:7]), axis=None),
                              theta=x[7])
    coords[idxs_to_shift] = rot_mat.dot(coords[idxs_to_shift].T).T

    # Compute the repulsive energy with the bonds crossing fragments
    energy = linker_energy(coords_i=coords[flat_x_motifs_idxs],
                           coords_j=coords[idxs_to_shift],
                           ij_bonds=ij_bonds)
    return energy


def x_motif_vec_and_fit(linker, fragment, template_linker, coords, x_motifs):
    """
    If this fragment contains an x motif then fit the coordinates will
    modify coords in place

    :param fragment:
    :param template_linker:
    :param coords:
    :param x_motifs:
    :return:
    """
    logger.info('Finding the X-motif shift vector for this fragment')

    for j, x_motif in enumerate(x_motifs):
        if all(atom_idx in fragment for atom_idx in x_motif.atom_ids):
            logger.info('Found fragment containing X-motif. Fitting coords')

            # Shifted template coordinates for this linker
            t_coords = np.array(template_linker.x_motifs[j].coords,
                                copy=True)

            x_coords = coords[np.array(x_motif.atom_ids, dtype=np.int)]
            # Centre on the average x, y, z coordinate
            x_centre = np.average(x_coords, axis=0)
            x_coords -= x_centre

            # Also centre the template coordinates
            t_centre = np.average(t_coords, axis=0)
            t_coords -= t_centre

            rot_mat = get_rot_mat_kabsch(p_matrix=x_coords,
                                         q_matrix=t_coords)

            # Apply the rotation to all the coords
            shift_x_coords = rot_mat.dot(x_coords.T).T

            # If the fitting cost of the x-motif fragment onto the template
            # is too high then try a dihedral rotation to reduce it
            if np.linalg.norm(shift_x_coords - t_coords) > 1:
                logger.warning('Poor fit. Needs rotation')
                rotate_x_motif_dihedral(linker, t_coords, x_motif, coords,
                                        fragment)
                # New X coordinates
                x_coords = coords[np.array(x_motif.atom_ids, dtype=np.int)]
                x_centre = np.average(x_coords, axis=0)
                x_coords -= x_centre

                rot_mat = get_rot_mat_kabsch(p_matrix=x_coords,
                                             q_matrix=t_coords)

            # Shift to the origin, rotate the fragment atoms then shift
            # back to the correct position
            f_idxs = np.array(fragment, dtype=np.int)

            coords[f_idxs] -= x_centre
            coords[f_idxs] = rot_mat.dot(coords[f_idxs].T).T
            coords[f_idxs] += t_centre

            return template_linker.x_motifs[j].norm_shift_vec

    return None


def fit_cost(t_coords, x_coords):
    """
    Calculate the fitting cost between two sets of coordinates

    :param t_coords: (np.ndarray)
    :param x_coords: (np.ndarray)
    :return:
    """

    # Centre on the average x, y, z coordinate
    x_centre = np.average(x_coords, axis=0)
    x_coords -= x_centre

    # Also centre the template coordinates
    t_centre = np.average(t_coords, axis=0)
    t_coords -= t_centre

    rot_mat = get_rot_mat_kabsch(p_matrix=x_coords,
                                 q_matrix=t_coords)

    # Apply the rotation to all the coords
    shift_x_coords = rot_mat.dot(x_coords.T).T

    return np.linalg.norm(shift_x_coords - t_coords)


def split_fragment_across_bond(linker, fragment, bond):
    """
    Split a set of atom indexes defining a fragment across a bond

    :param linker:
    :param fragment:
    :param bond:
    :return:
    """

    graph = linker.graph.copy()
    graph.remove_edge(*bond)
    components = list(nx.connected_components(graph))

    if len(components) != 2:
        raise ValueError('Failed to split across bond - likely in a ring')

    # Atom indexes in the left and right portion of the fragment
    return tuple([i for i in fragment if i in components[j]] for j in range(2))


def rotate_x_motif_dihedral(linker, template_coords, x_motif,
                            coords, fragment):
    """
    Maximise the fit between linker X-motif with some coordinates and the
    template coordinates by rotating around a dihedral contained within the
    X-motif. Only rotate on single bonds

    :param linker: (cgbind.linker.Linker)
    :param template_coords:
    :param x_motif:
    :param coords:
    :param fragment: (list(int)) Atom indexes that need to be modified
    """
    logger.info('Attempting dihedral X-motif rotation to improve the fit')

    try:
        bonds = linker.get_single_bonds()
    except ValueError:
        logger.warning('Could not rotate - single bonds not defined')
        return None

    for (i, j) in bonds:
        if i in x_motif.atom_ids and j in x_motif.atom_ids:
            logger.info(f'Found rotatable bond in X-motif {i}-{j}')

            try:
                l_idxs, r_idxs = split_fragment_across_bond(linker, fragment,
                                                            bond=(i, j))
            except ValueError:
                logger.warning('Could not split and rotate')
                return None

            origin = np.array(coords[i], copy=True)
            axis = coords[i] - coords[j]

            min_cost, best_rot_mat, best_idxs = None, None, None

            for rot_idxs in (l_idxs, r_idxs):
                # Convert to an array that can be used for indexing
                rot_idxs = np.array(rot_idxs, dtype=np.int)

                x_idxs = np.array(x_motif.atom_ids, dtype=np.int)

                # Possible rotations for single bonds
                for theta in [0, 120, 180, 240]:
                    rot_coords = np.array(coords, copy=True)
                    rot_mat = rotation_matrix(axis, theta=np.deg2rad(theta))

                    # Shift to the origin, rotate, and shift back
                    rot_coords -= origin
                    rot_coords[rot_idxs] = rot_mat.dot(rot_coords[rot_idxs].T).T

                    cost = fit_cost(t_coords=template_coords,
                                    x_coords=rot_coords[x_idxs])

                    if min_cost is None or cost < min_cost:
                        min_cost = cost
                        best_rot_mat = rot_mat
                        best_idxs = rot_idxs

            logger.info(f'Best indexes to rotate: {best_idxs}')
            logger.info(f'Lowest fitting cost = {min_cost:.4f}')

            # Apply the rotation for real
            coords[best_idxs] -= origin
            coords[best_idxs] = best_rot_mat.dot(coords[best_idxs].T).T
            coords[best_idxs] += origin

    return None
