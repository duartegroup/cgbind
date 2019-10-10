from cgbind.log import logger
from copy import deepcopy
import numpy as np
from cgbind.geom import rotation_matrix
from cgbind.geom import xyz2coord
from cgbind.geom import calc_com
from cgbind.geom import calc_normalised_vector
from cgbind.geom import cat_cage_subst_coords


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

    # Construct a set of orthogonal vectors with the first
    rot1 = calc_normalised_vector(np.zeros(3), cage_coords[cage.m_ids[0]])
    rot2 = calc_normalised_vector(rot1[0] * rot1, np.array([1.0, 0.0, 0.0]))
    rot3 = calc_normalised_vector(rot1[1] * rot1 + rot2[1] * rot2, np.array([0.0, 1.0, 0.0]))

    subst_coords = rot_minimise_repulsion(subst_coords, cage_coords, rot_axes=[rot1, rot2, rot3], n_rot_steps=100)

    xyzs = cat_cage_subst_coords(cage, substrate, cage_coords, subst_coords)

    return xyzs


def rot_minimise_repulsion(subst_coords, cage_coords, rot_axes, n_rot_steps=100):
    """
    Rotate a substrate around the M-M (z) axis as to reduce the steric repulsion between it and the cage
    :param subst_coords:
    :param cage_coords:
    :param rot_axes: List of orthogonal normalised vectors
    :param n_rot_steps: Number of rotation steps to perform in each axis
    :return:
    """
    logger.info('Minimising steric repulsion between substrate and cage by rotation')

    best_theta = 0.0
    best_sum_cage_substrate_inverse_dists = 9999999.9
    best_axis = rot_axes[0]                         # TODO allow for a combination of axes

    for rot_axis in rot_axes:
        for theta in np.linspace(0, np.pi / 2.0, n_rot_steps):
            tmp_substrate_coords = deepcopy(subst_coords)
            rot_matrix = rotation_matrix(rot_axis, theta)
            tmp_rot_substrate_coords = [np.matmul(rot_matrix, coord) for coord in tmp_substrate_coords]

            sum_cage_substrate_inverse_dists = 0.0
            for cage_coord in cage_coords:
                for substrate_coord in tmp_rot_substrate_coords:
                    dist = np.linalg.norm(cage_coord - substrate_coord)
                    if dist < 2.0:
                        sum_cage_substrate_inverse_dists += 100.0 / dist  # Penalise very short distances

            if sum_cage_substrate_inverse_dists < best_sum_cage_substrate_inverse_dists:
                best_theta = theta
                best_sum_cage_substrate_inverse_dists = sum_cage_substrate_inverse_dists
                best_axis = rot_axis

    rot_matrix = rotation_matrix(best_axis, best_theta)
    rot_substrate_coords = [np.matmul(rot_matrix, coord) for coord in subst_coords]

    return rot_substrate_coords

