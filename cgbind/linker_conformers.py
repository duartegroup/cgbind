import itertools
import networkx as nx
import numpy as np
from cgbind.exceptions import CgbindCritical
from cgbind.geom import rotation_matrix
from cgbind.log import logger


def set_best_fit_linker(linker, x_motifs, template_linker):
    """
    Generate conformers of a linker in the space of single bond dihedrals that
    alter the fit to the template. This should restrict the space considerably
    enabling the enumeration of all possibilities

    :param linker: (cgbind.linker.Linkler)
    :param x_motifs: (list(cgbind.x_motifs.Xmotif))
    :param template_linker: (cgbind.templates.Linker)
    """
    logger.info('Generating conformers of the linker that change the template '
                'fit')

    t_coords = np.vstack(tuple(xm.coords for xm in template_linker.x_motifs))

    # Dictionary of bonds as the tuple of atom indexes and the indexes of atoms
    # to rotate on the left and right of the split
    bonds_and_lr_idxs = get_rot_bonds_and_idxs(linker, t_coords, x_motifs)

    # Possible angles for the rotatable dihedrals
    thetas_list = itertools.product([0, 2*np.pi/3, np.pi, 4*np.pi/3],
                                    repeat=len(bonds_and_lr_idxs))

    # Indexes in the coordinates that will be fit to the template
    fit_idxs = [atom_id for motif in x_motifs for atom_id in motif.atom_ids]

    # Calculate the best list of thetas using the Cython extension
    try:
        from conf_fit import opt_coords
    except ModuleNotFoundError:
        return CgbindCritical('opt_coords module not built')

    best_coords = opt_coords(np.array(list(thetas_list)),
                             bonds_and_lr_idxs,
                             py_coords=linker.get_coords(),
                             py_template_coords=t_coords,
                             fit_idxs=np.array(fit_idxs, dtype=np.int))

    linker.set_atoms(coords=best_coords)
    return None


def get_rot_bonds_and_idxs(linker, t_coords, x_motifs):
    """
    Find the bonds in the linker that can be rotated to alter the fitting cost

    :param linker:
    :param t_coords: (np.ndarray) Template coordinates
    :param x_motifs:
    :return:
    """
    bonds = linker.get_single_bonds()
    graph = linker.graph
    coords = linker.get_coords()

    bonds_and_rot_idxs = {}

    for (i, j) in bonds:

        graph.remove_edge(i, j)
        components = list(nx.connected_components(graph))
        graph.add_edge(i, j)

        # Only consider dihedrals that aren't in rings
        if len(components) != 2:
            continue

        l_idxs, r_idxs = components
        # Only consider dihedrals between large-ish fragments
        if len(l_idxs) < 4 or len(r_idxs) < 4:
            continue

        # Only consider dihedrals that upon rotation lead to a change in the
        # fitting cost upon rotation
        rot_idxs = np.array(list(l_idxs)
                            , dtype=np.int)
        if delta_cost_on_rot(coords, t_coords, x_motifs, rot_idxs,
                             bond=(i, j), theta=np.pi/2) < 1E-3:
            continue

        # TODO: Is ok to only consider one side??
        bonds_and_rot_idxs[(i, j)] = l_idxs

    logger.info(f'Have *{len(bonds_and_rot_idxs)}* dihedrals that can be '
                f'rotated')
    return bonds_and_rot_idxs


def delta_cost_on_rot(coords, t_coords, x_motifs, rot_idxs, bond, theta):
    """
    Rotate some coordinates along a bond axis an amount theta radians and
    return the difference in fit to the  between

    :param coords:
    :param t_coords:
    :param x_motifs:
    :param rot_idxs:
    :param bond:
    :param theta:
    :return:
    """
    logger.info('Calculating ∆C for the fitting cost (C) when this bond'
                ' is rotated')

    x_idxs = np.array([atom_id for motif in x_motifs
                       for atom_id in motif.atom_ids], dtype=np.int)
    curr_cost = fit_cost(coords[x_idxs], t_coords)

    rot_coords = coords_rotated_dihedral(coords, bond, theta, rot_idxs,
                                         copy=True)

    rot_cost = fit_cost(rot_coords[x_idxs], t_coords)
    delta_cost = np.abs(curr_cost - rot_cost)

    logger.info(f'∆C = {delta_cost:.4f}')
    return delta_cost


def coords_rotated_dihedral(coords, bond, theta, rot_idxs, copy=True):
    """
    Generate coordinates by rotation around a bond

    :param coords:
    :param bond:
    :param theta: (float)
    :param rot_idxs:
    :param copy:
    :return:
    """
    rot_coords = np.array(coords, copy=copy)

    l_idx, r_idx = bond
    rot_mat = rotation_matrix(axis=coords[l_idx] - coords[r_idx],
                              theta=theta)

    # Rotate the portion of the structure to rotate about the bonding axis
    rot_coords[rot_idxs] -= coords[l_idx]
    rot_coords[rot_idxs] = rot_mat.dot(rot_coords[rot_idxs].T).T
    rot_coords[rot_idxs] += coords[l_idx]

    return rot_coords


def fit_cost(x_coords, t_coords):
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

    rot_mat = kabsch_rot_matrix(p_matrix=x_coords,
                                q_matrix=t_coords)

    # Apply the rotation to all the coords
    shift_x_coords = rot_mat.dot(x_coords.T).T

    return np.linalg.norm(shift_x_coords - t_coords)


def kabsch_rot_matrix(p_matrix, q_matrix):
    """
    Get the optimal rotation matrix with the Kabsch algorithm. Notation is from
    https://en.wikipedia.org/wiki/Kabsch_algorithm

    :param p_matrix: (np.ndarray)
    :param q_matrix: (np.ndarray)
    :return: (np.ndarray) rotation matrix
    """

    h = np.matmul(p_matrix.transpose(), q_matrix)
    u, s, vh = np.linalg.svd(h)
    d = np.linalg.det(np.matmul(vh.transpose(), u.transpose()))
    int_mat = np.identity(3)
    int_mat[2, 2] = d
    rot_matrix = np.matmul(np.matmul(vh.transpose(), int_mat), u.transpose())

    return rot_matrix