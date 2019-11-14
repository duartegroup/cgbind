from cgbind.log import logger
from copy import deepcopy
import numpy as np
from scipy.optimize import minimize, Bounds
from scipy.spatial import distance_matrix
from cgbind import geom
from cgbind.atoms import heteroatoms
from cgbind.geom import rotation_matrix
from cgbind.geom import xyz2coord
from cgbind.geom import calc_com


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


def add_substrate_com(cage, substrate):
    """
    Add a substrate the centre of a cage defined by its centre of mass (com)
    :param cage: A Cage object
    :param substrate: A Substrate object
    :return:
    """
    logger.info('Adding substrate to the cage COM and minimising repulsion')

    cage_coords = xyz2coord(cage.xyzs)
    subst_coords = xyz2coord(substrate.xyzs)

    subst_com = calc_com(substrate.xyzs)

    # Get the centroid of the metals, which will is ~ the cage COM for a symmetric system
    metal_coords = np.array([cage_coords[i] for i in cage.m_ids])
    centroid = np.average(metal_coords, axis=0)

    # Shift both the cage and the substrate to the origin
    cage_coords = [coord - centroid for coord in cage_coords]
    subst_coords = [coord - subst_com for coord in subst_coords]
    subst_x_ids = [i for i in range(substrate.n_atoms) if substrate.xyzs[i][0] in heteroatoms]

    logger.info('Minimising steric repulsion between substrate and cage & metal-X stom dists  by rotation')
    min_result = minimize(cost_repulsion_and_x_metal_dist, x0=np.array([1.0, 1.0, 1.0]),
                          args=(subst_coords, cage_coords, cage.m_ids, subst_x_ids), method='L-BFGS-B',
                          bounds=Bounds(lb=0.0, ub=2*np.pi))

    logger.info(f'Optimum rotation is {min_result.x} rad in x, y, z')
    subst_coords = cost_repulsion_and_x_metal_dist(min_result.x, subst_coords, cage_coords,
                                                   cage.m_ids, subst_x_ids, return_cost=False)

    xyzs = cat_cage_subst_coords(cage, substrate, cage_coords, subst_coords)

    return xyzs


def cost_repulsion_and_x_metal_dist(x, subst_coords, cage_coords, metal_ids, subst_x_ids, return_cost=True):
    """
    Calculate the cost function for a particular x, which contains the rotations in x, y, z cartesian directions
    by rotating the subst_coords with sequencial rotations. The cost function is 1/r^6 between the cage and substrate
    and tanh(r - 5.0) + 1.0 for the substrate heteroatoms – metal atoms.

    :param x: (np.ndarray)
    :param subst_coords: (list(np.ndarray))
    :param cage_coords: (list(np.ndarray))
    :param metal_ids: (list(int))
    :param subst_x_ids: (list(int))
    :param return_cost: (bool)
    """

    x_rot, y_rot, z_rot = x

    rot_matrix = np.identity(3)
    rot_matrix = np.matmul(rot_matrix, rotation_matrix(axis=geom.i, theta=x_rot))
    rot_matrix = np.matmul(rot_matrix, rotation_matrix(axis=geom.j, theta=y_rot))
    rot_matrix = np.matmul(rot_matrix, rotation_matrix(axis=geom.k, theta=z_rot))

    rot_substrate_coords = [np.matmul(rot_matrix, coord) for coord in deepcopy(subst_coords)]

    inv_distance_mat = np.power(distance_matrix(rot_substrate_coords, cage_coords), -6)

    # Compute the distances from the metals to the substrate heteroatoms
    metal_subt_x_dists = np.array([np.linalg.norm(rot_substrate_coords[i] - cage_coords[j])
                                   for i in subst_x_ids for j in metal_ids])

    # Weight the distances so short distances count  little to the cost function
    weighted_subt_x_dists = np.tanh(metal_subt_x_dists - 5.0) + 1.0

    cost = np.sum(inv_distance_mat) + np.sum(weighted_subt_x_dists)
    if return_cost:
        return cost

    return rot_substrate_coords

