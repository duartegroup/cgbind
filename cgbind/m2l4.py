import numpy as np
from copy import deepcopy
from .log import logger
from .input_output import print_output
from .geom import rotation_matrix
from .geom import xyz2coord
from .geom import calc_midpoint
from .geom import get_closest_bonded_atom_id
from .geom import calc_normalised_vector
from .geom import calc_com
from .geom import calc_dist
from .geom import cat_cage_subst_coords


def get_best_linker_conformer(linker, cavity_size=2.2):
    """
    For a linker to be suitable for use in a Pd2l4 cage of the sort
            N
           /|
    --N - M- N--
    |    /  |  |
    |   N   |  |
    |  |    |  |
    |  |    N  |
    |  |   /   |
    --N|- M- N--
       |/
       N

    the confomer of the linker MUST have a cavity between the midpoint of the two atoms. This function will process
    the conformers and return the ids of those which meet the criteria. The N atoms coordinated to the Pd are found
    by searching for the longest N-N distance. The most suitable (best) confomer will be that whose N–N distance
    is minimised, which will occur if they lie in the same plane.

    :param linker: A linker object
    :param cavity_size: Radius of the sphere surrounding the midpoint of the N-N vector. If any atoms are within
    this distance the confomer is not suitable.
    :return: Best confomer id
    """
    logger.info('Getting the best linker conformer for {}'.format(linker.name))
    logger.warning('Will only consider N atoms as donor atoms')  # TODO extend to other possible donor atoms

    if linker.conf_ids is None or linker.conf_xyzs is None:
        return None

    best_conf_id = None
    best_conf_n_n_dist = 9999.9

    for conf_id in linker.conf_ids:

        xyzs = linker.conf_xyzs[conf_id]
        n1_atom_id, n2_atom_id = get_furthest_away_nitrogen_atoms(xyzs)
        if n1_atom_id is None or n2_atom_id is None:
            logger.error('Linker did not have two nitrogen atoms')
            return None

        n1_n2_dist = calc_dist(xyzs, atom_i=n1_atom_id, atom_j=n2_atom_id)

        x_x_midpoint = (xyz2coord(xyzs[n1_atom_id]) + xyz2coord(xyzs[n2_atom_id])) / 2.0
        suitable_confomer = True

        for line_i in xyzs:
            dist_atom_midpoint = np.linalg.norm(xyz2coord(line_i) - x_x_midpoint)
            if dist_atom_midpoint < cavity_size:
                suitable_confomer = False

        if suitable_confomer:
            if n1_n2_dist < best_conf_n_n_dist:
                best_conf_id = conf_id
                best_conf_n_n_dist = n1_n2_dist
                linker.x_atom_ids = [n1_atom_id, n2_atom_id]

    if best_conf_id is None:
        logger.warning('Could not find a best conformer for {}'.format(linker.name))
        return None

    print_output('Best conformer of', linker.name, 'Found')
    logger.info('Found a suitable conformer with id {}'.format(best_conf_id))
    return best_conf_id


def get_furthest_away_nitrogen_atoms(xyzs):
    """
    Get the pair of nitrogen atoms that are the furthest away in space. This is perhaps not the best way to find the
    donor nitrogen atoms...

    :param xyzs: (list(list))
    :return: (tuple) atom indexes
    """

    max_n_n_dist = 0.0
    n1_atom_id, n2_atom_id = None, None

    for i, line_i in enumerate(xyzs):
        for j, line_j in enumerate(xyzs):
            if i > j and line_i[0] == 'N' and line_j[0] == 'N':
                dist = np.linalg.norm(xyz2coord(line_i) - xyz2coord(line_j))
                if dist > max_n_n_dist:
                    max_n_n_dist = dist
                    n1_atom_id, n2_atom_id = i, j

    return n1_atom_id, n2_atom_id


def get_mid_atom_id(linker):
    """
    For a linker then the 'mid_atom' is needed to define a plane from which the M2L4 cage may be built. This is
    found by calculating the abs difference in the A-X1, A-X2 distances for all atoms A, the minimum of which will
    be the atom closest to the midpoint of the linker as possible.
    :return: The ID of the atom in the xyz list
    """
    logger.info('Calculating mid atom id for linker {}'.format(linker.name))

    if linker.xyzs is None:
        logger.error('Linker xyzs were None. Cannot find the mid atom id')
    if linker.x_atom_ids is None:
        logger.error('Could not find the mid atom. Have no x_atom_ids')
        return None

    mid_atom_id = 0
    x1_atom_id, x2_atom_id = linker.x_atom_ids
    x1_coords, x2_coords = xyz2coord(linker.xyzs[x1_atom_id]), xyz2coord(linker.xyzs[x2_atom_id])
    linker.x_x_midpoint = calc_midpoint(x1_coords, x2_coords)

    min_dist_diff = 999.9

    for i in range(len(linker.xyzs)):
        atom_i_coords = xyz2coord(linker.xyzs[i])
        dist_diff = np.abs(np.linalg.norm(atom_i_coords - x1_coords) - np.linalg.norm(atom_i_coords - x2_coords))
        if dist_diff < min_dist_diff:
            min_dist_diff = dist_diff
            mid_atom_id = i

    return mid_atom_id


def build(cage, linker):
    """
    From a linker object rotate the linker around a defined center to create a L4 scaffold, then add the M atoms
    to the midpoints between the N-atoms from 'opposite' linkers
    :param cage: A Cage object
    :param linker: A linker object
    :return: The xyzs of the full cage
    """
    logger.info('Building M2L4 metallocage'.format(cage.name))

    if any(k is None for k in (linker.mid_atom_id, linker.xyzs, linker.x_x_midpoint, linker.x_atom_ids)):
        return None

    x_x_dist = 4.1  # X - X distance in a N-M-N unit
    delta_x_x_dist = 999.9
    xyzs = []
    trns_dist = -20.0

    while delta_x_x_dist > 0.05:

        xyzs = []

        # ------ Add the 4 linkers in the correct orientation ----------

        trns_dist += 0.05
        cage.centre = get_cage_centre(linker, trns_dist)

        first_linker_coords = xyz2coord(linker.xyzs) - cage.centre

        midpoint_n1_vector = linker.x_x_midpoint - xyz2coord(linker.xyzs[linker.x_atom_ids[0]])
        norm_midpoint_n1_vector = midpoint_n1_vector / np.linalg.norm(midpoint_n1_vector)

        rotation_axis = norm_midpoint_n1_vector

        for angle in [0.0, np.pi / 2.0, np.pi, 3.0 * np.pi / 2.0]:
            linker_rot_xyzs = []
            for i in range(len(linker.xyzs)):
                linker_coord = np.matmul(rotation_matrix(rotation_axis, angle), first_linker_coords[i])
                linker_xyz = [linker.xyzs[i][0]] + linker_coord.tolist()
                linker_rot_xyzs.append(linker_xyz)

            xyzs += linker_rot_xyzs

        actual_n_n_dist = np.linalg.norm(xyz2coord(xyzs[linker.x_atom_ids[0]]) -
                                         xyz2coord(xyzs[linker.x_atom_ids[0] + 2 * linker.n_atoms]))

        def check_linker_geom(tolerance=0.8):
            """
            All atoms in opposite linkers need to be further away than the X - X distance
            :return: A large val which will be added to delta_x_x_dist, ensuring the while loop isn't exited
            """

            for i in range(len(linker.xyzs)):
                for j in range(len(linker.xyzs)):
                    if i > j:
                        dist = np.linalg.norm(xyz2coord(xyzs[i]) - xyz2coord(xyzs[j + 2 * linker.n_atoms]))
                        if dist + tolerance < actual_n_n_dist:
                            return 99.9

            return 0.0

        delta_x_x_dist = np.abs(actual_n_n_dist - x_x_dist) + check_linker_geom()

        if trns_dist > 100:
            cage.reasonable_geometry = False
            logger.warning("No sensible geometry can be found for {}".format(linker.name))
            return None

    # ------ Add the M atoms ----------

    for i in range(2):
        x_x_midpoint_oposite_linkers = calc_midpoint(xyz2coord(xyzs[linker.x_atom_ids[i]]),
                                                     xyz2coord(xyzs[linker.x_atom_ids[i] + 2 * linker.n_atoms]))
        # Add the metal atom to the xyz list. Located at the midpoint between 2 N atoms on opposite sides
        xyzs.append([cage.metal] + x_x_midpoint_oposite_linkers.tolist())

    return xyzs


def get_cage_centre(linker, trns=3.0):
    """
    Get the centre defined as the following

           -------X
           |
           |
    midatom|             centre
           |
           |
           -------X
    :param linker: A Linker object
    :param trns: Distance from the mid atom to the centre
    :return: Coordinates of the centre
    """

    mid_atom_coords = xyz2coord(linker.xyzs[linker.mid_atom_id])
    midpoint_mid_atom_vector = linker.x_x_midpoint - mid_atom_coords
    norm_midpoint_mid_atom_vector = mid_atom_coords / np.linalg.norm(midpoint_mid_atom_vector)

    return linker.x_x_midpoint - trns * norm_midpoint_mid_atom_vector


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


def add_substrate_x_x(cage, substrate):
    """
    Add the substrate xyzs to the list of cage xyzs. Locate so the X-X vector is aligned with the M-M vector,
    then rotate in the M-M axis to minimise steric clash between the substrate and the cage

              N
             /|
    ----N - M- N----
    |      /  |    |
    |     N X |    |
    |   s|ubstate  |
    |    |  X N    |
    |    |   /     |
    ----N|- M- N----
         |/
         N

    :param cage: A Cage object
    :param substrate: A Substrate object
    :return: List of xyzs with both cage and substrate xyzs
    """
    logger.info('Adding substrate with 2(?) heteroatoms aligned along the M–M vector')

    cage_coords = xyz2coord(cage.xyzs)
    subst_coords = xyz2coord(substrate.xyzs)
    subst_x_id_0, subst_x_id_1 = substrate.x_x_atom_ids[0], substrate.x_x_atom_ids[1]
    substrate_center = calc_midpoint(subst_coords[subst_x_id_1], subst_coords[subst_x_id_0])

    # Translate the substrate to the M-M midpoints

    m_m_midpoint = calc_midpoint(cage_coords[cage.m_ids[1]], cage_coords[cage.m_ids[0]])
    cage_coords = [coord - m_m_midpoint for coord in cage_coords]  # Centre the cage on the M-M midpoint
    subst_coords = [coord - (m_m_midpoint + substrate_center) for coord in subst_coords]

    # Rotate the substrate so the X–X vector of the substrate is aligned with the M-M vector

    substrate_centre_x_vector = calc_normalised_vector(m_m_midpoint, subst_coords[subst_x_id_0])
    cage_m_center_vector = calc_normalised_vector(m_m_midpoint, cage_coords[cage.m_ids[1]])

    theta = np.arccos(np.dot(substrate_centre_x_vector, cage_m_center_vector))
    normal = np.cross(substrate_centre_x_vector, cage_m_center_vector)
    rot_matrix = rotation_matrix(normal, theta)
    subst_coords = [np.matmul(rot_matrix, coord) for coord in subst_coords]

    rot_substrate_coords = rot_minimise_repulsion(subst_coords, cage_coords, rot_axes=[cage_m_center_vector])
    xyzs = cat_cage_subst_coords(cage, substrate, cage_coords, rot_substrate_coords)

    return xyzs


def add_substrate_x(cage, substrate, m_x_dist=3.5):
    """
    Add a substrate with at least one heteroatom (X) which can form C(2)H hydrogen bonds, pointing at a M atom
    will just pick the first reasonable conformer...

    :param cage: A cage object
    :param substrate: A substrate object
    :param m_x_dist: Guess at the M-X bond distance in Å
    :return:
    """
    logger.info('Adding substrate with 1(?) heteroatom aligned along the M–M vector')

    cage_coords = xyz2coord(cage.xyzs)
    subst_coords = xyz2coord(substrate.xyzs)

    for x_atom_id in substrate.x_atom_ids:

        # Shift the M and X atoms to (0.0, 0.0, 0.0), shift the X atom to the first M atom
        cage_coords = [coord - cage_coords[cage.m_ids[0]] for coord in cage_coords]
        subst_coords = [coord - subst_coords[x_atom_id] for coord in subst_coords]

        # Assume the X atom is monovalent.... and bonded to an element E
        e_atom_id = get_closest_bonded_atom_id(substrate.xyzs, x_atom_id)
        e_x_bond_vector = calc_normalised_vector(subst_coords[e_atom_id], subst_coords[x_atom_id])
        m_m_vector = calc_normalised_vector(cage_coords[cage.m_ids[1]], cage_coords[cage.m_ids[0]])

        # Rotate the substrate to the E-X bond is inline with the M-M vector
        theta = np.arccos(np.dot(e_x_bond_vector, m_m_vector))
        normal = np.cross(e_x_bond_vector, m_m_vector)
        rot_matrix = rotation_matrix(normal, theta)
        subst_coords = [np.matmul(rot_matrix, coord) for coord in subst_coords]

        # Translate the substrate so the M-X distance is as defined
        subst_coords = [coord - (m_x_dist * m_m_vector) for coord in subst_coords]

        # Rotate in the M-M axis to minimise the steric repulsion
        subst_coords = rot_minimise_repulsion(subst_coords, cage_coords, rot_axes=[m_m_vector])

        xyzs = cat_cage_subst_coords(cage, substrate, cage_coords, subst_coords)

        # TODO loop over possible binding atoms, perhaps printing all possible binding modes

        return xyzs


def add_substrate_com(cage, substrate):
    """
    Add a substrate the centre of a cage defined by its centre of mass (com)
    :param cage: A Cage object
    :param substrate: A Substrate object
    :return:
    """
    logger.info('Adding substrate with 0(?) heteroatoms to the cage COM')

    cage_coords = xyz2coord(cage.xyzs)
    subst_coords = xyz2coord(substrate.xyzs)

    subst_com = calc_com(substrate.xyzs)

    m_m_midpoint = calc_midpoint(cage_coords[cage.m_ids[1]], cage_coords[cage.m_ids[0]])
    subst_coords = [coord - (m_m_midpoint + subst_com) for coord in subst_coords]

    rot1 = calc_normalised_vector(m_m_midpoint + subst_com, cage_coords[cage.m_ids[0]])
    rot2 = calc_normalised_vector(rot1[0] * rot1, np.array([1.0, 0.0, 0.0]))
    rot3 = calc_normalised_vector(rot1[1] * rot1 + rot2[1] * rot2, np.array([0.0, 1.0, 0.0]))

    subst_coords = rot_minimise_repulsion(subst_coords, cage_coords, rot_axes=[rot1, rot2, rot3], n_rot_steps=100)

    xyzs = cat_cage_subst_coords(cage, substrate, cage_coords, subst_coords)

    return xyzs
