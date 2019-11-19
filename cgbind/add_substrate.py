from cgbind.log import logger
from copy import deepcopy
import numpy as np
from scipy.optimize import minimize, Bounds
from scipy.spatial import distance_matrix
from cgbind import geom
from cgbind.atoms import get_vdw_radii
from cgbind.geom import rotation_matrix
from cgbind.geom import xyz2coord
from cgbind.geom import calc_com


def cage_subst_repulsion_func(cage, substrate, cage_coords, subst_coords):

    dist_mat = distance_matrix(cage_coords, subst_coords)

    exponents = np.add.outer(np.array(substrate.vdw_radii),
                             np.array(cage.vdw_radii))

    # TODO work out this function
    energy_mat = np.exp((1 - dist_mat / exponents))
    energy = np.sum(energy_mat)

    #       Energy for a binding cage will should be negative so this magic number controls if it is
    return energy - 2


def add_substrate_com(cagesubt):
    """
    Add a substrate the centre of a cage defined by its centre of mass (com)
    :param cagesubt: (object)
    :return:
    """
    logger.info('Adding substrate to the cage COM and minimising the energy')

    # Minimum energy initialisation and the x parameter array (angles to rotate about the x, y, z axes)
    min_energy, curr_x = 9999999999.9, np.zeros(3)

    # Optimum (minimum energy) conformer
    best_i = 0

    c, s = cagesubt.cage, cagesubt.substrate
    cage_coords = get_centered_cage_coords(c.xyzs, c.m_ids)
    c.vdw_radii = [get_vdw_radii(atom_label=xyz[0]) for xyz in c.xyzs]

    if cagesubt.n_subst_confs > 1:
        s.gen_confs(n_confs=cagesubt.n_subst_confs)

    for i, substrate_xyzs in enumerate(s.conf_xyzs):
        subst_coords = get_centered_substrate_coords(substrate_xyzs)
        s.vdw_radii = [get_vdw_radii(atom_label=xyz[0]) for xyz in s.xyzs]

        for _ in range(cagesubt.n_init_geom):
            rot_angles = 2.0 * np.pi * np.random.rand(3)        # rand generates in [0, 1] so multiply with

            # Minimise the energy with a BFGS minimiser supporting bounds on the values (rotation is periodic)
            result = minimize(get_energy, x0=np.array(rot_angles),
                              args=(c, s, cagesubt.energy_func, subst_coords, cage_coords), method='L-BFGS-B',
                              bounds=Bounds(lb=0.0, ub=2*np.pi))

            energy = result.fun
            logger.info(f'Energy is {energy}')

            if energy < min_energy:
                min_energy = energy
                curr_x = result.x
                best_i = i

    # Get the new coordinates of the substrate that has been rotated appropriately
    new_subst_coords = get_rotated_subst_coords(curr_x, subst_coords=xyz2coord(s.conf_xyzs[best_i]))

    return cat_cage_subst_coords(c, s, cage_coords, new_subst_coords)


def get_centered_cage_coords(cage_xyzs, cage_m_ids):
    """Get the cage coordinates that had been translated to the cage centroid"""

    cage_coords = xyz2coord(cage_xyzs)
    metal_coords = np.array([cage_coords[i] for i in cage_m_ids])
    centroid = np.average(metal_coords, axis=0)

    return [coord - centroid for coord in cage_coords]


def get_centered_substrate_coords(substrate_xyzs):
    """Get the substrate coorindated that have been translated to its center of mass"""

    subst_coords = xyz2coord(substrate_xyzs)
    subst_com = calc_com(substrate_xyzs)
    return [coord - subst_com for coord in subst_coords]


def cat_cage_subst_coords(cage, substrate, cage_coords, substrate_coords):
    """
    Catenate some coordinates into a set of xyzs by adding back the atom labels from the original xyzs

    :param cage:
    :param substrate:
    :param cage_coords:
    :param substrate_coords:
    :return:
    """
    logger.info('Appending substrate coordinates to cage coordinates')

    xyzs = [[cage.xyzs[n][0]] + cage_coords[n].tolist() for n in range(len(cage.xyzs))]
    cage.substrate_atom_ids = list(range(len(xyzs), len(xyzs) + len(substrate.xyzs)))
    xyzs += [[substrate.xyzs[n][0]] + substrate_coords[n].tolist() for n in range(len(substrate.xyzs))]

    return xyzs


def get_rotated_subst_coords(x, subst_coords):
    """Get substrate coordinates that have been rotated by x[0] rad in the x axis etc."""

    x_rot, y_rot, z_rot = x

    rot_matrix = np.identity(3)
    rot_matrix = np.matmul(rot_matrix, rotation_matrix(axis=geom.i, theta=x_rot))
    rot_matrix = np.matmul(rot_matrix, rotation_matrix(axis=geom.j, theta=y_rot))
    rot_matrix = np.matmul(rot_matrix, rotation_matrix(axis=geom.k, theta=z_rot))

    return [np.matmul(rot_matrix, coord) for coord in deepcopy(subst_coords)]


def get_energy(x, cage, substrate, energy_func, cage_coords, subst_coords):
    """
    Calculate the cost function for a particular x, which contains the rotations in x, y, z cartesian directions
    by rotating the subst_coords with sequencial rotations. The cost function is 1/r^6 between the cage and substrate
    and tanh(r - 5.0) + 1.0 for the substrate heteroatoms – metal atoms.
    """
    rot_substrate_coords = get_rotated_subst_coords(x, subst_coords)
    energy = energy_func(cage, substrate, cage_coords, rot_substrate_coords)

    return energy


cage_subst_repulsion_func.__name__ = 'repulsion'
energy_funcs = [cage_subst_repulsion_func]
