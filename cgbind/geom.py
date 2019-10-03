import numpy as np
from cgbind.log import logger
from cgbind.atoms import avg_bond_lengths
from cgbind.atoms import get_atomic_mass


def calc_com(xyzs):
    """
    Calculate the centre of mass for a list of xyzs
    :param xyzs:
    :return:
    """
    logger.info('Calculating centre of mass ')

    com = np.zeros(3)
    total_mass = 0.0

    for n in range(len(xyzs)):
        atom_mass = get_atomic_mass(atom_label=xyzs[n][0])
        total_mass += atom_mass

        com += atom_mass * xyz2coord(xyzs[n])

    return com / total_mass


def calc_cdist(coord1, coord2):
    return np.linalg.norm(coord1 - coord2)


def calc_normalised_vector(coord1, coord2):

    vec = coord2 - coord1
    return vec / np.linalg.norm(vec)


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


def get_closest_bonded_atom_id(xyzs, atom_id, tolerance=0.1, default_bond_length=1.5):
    """
    Get the closest atom to a particular id
    :param xyzs:
    :param atom_id:
    :param tolerance: Tolerance on deviation from avg bond length. Default 10%
    :param default_bond_length:
    :return:
    """
    logger.info('Getting closest bonded atom id to {}'.format(atom_id))

    min_dist = 999.9
    closest_bonded_atom_id = 0

    for n in range(len(xyzs)):
        if n != atom_id:
            dist = np.linalg.norm(xyz2coord(xyzs[n]) - xyz2coord(xyzs[atom_id]))
            if dist < min_dist:
                atoms_key1, atoms_key2 = (xyzs[n][0] + xyzs[atom_id][0]), (xyzs[atom_id][0] + xyzs[n][0])
                if atoms_key1 in avg_bond_lengths.keys():
                    max_bond_length = (1.0 + tolerance) * avg_bond_lengths[atoms_key1]
                elif atoms_key2 in avg_bond_lengths.keys():
                    max_bond_length = (1.0 + tolerance) * avg_bond_lengths[atoms_key2]
                else:
                    logger.warning('Couldn\'t find {}–{} avg. bond distance in dict'.format(atoms_key1, atoms_key2))
                    max_bond_length = default_bond_length

                if dist < max_bond_length:
                    min_dist = dist
                    closest_bonded_atom_id = n

    return closest_bonded_atom_id


def get_ids_max_dist(distance_matrix):
    """
    Poor mans np.argmax
    :param distance_matrix:
    :return:
    """

    max_dist = 0.0
    atom_ids = [0, 0]
    for atom_i in range(len(distance_matrix)):
        for atom_j in range(len(distance_matrix)):
            if distance_matrix[atom_i, atom_j] > max_dist:
                atom_ids = [atom_i, atom_j]
                max_dist = distance_matrix[atom_i, atom_j]

    return atom_ids


def get_neighbour(atom_id, distance_matrix, proximity):
    """
    Get the  neighbour to a specific atom given a distance matrix
    :param atom_id:
    :param distance_matrix:
    :param proximity: Proximity to the atom_id i.e. 1 = nearest neighbour, 2 = next-nearest neighbour
    :return: Id of the nearest neighbour to an atom
    """
    # logger.info('Getting {}th nearest neighbour to atom id {}'.format(proximity, atom_id))
    neighbour = 0
    dists_low_to_high = sorted(distance_matrix[atom_id])

    for atom_j in range(len(distance_matrix[atom_id])):
        if distance_matrix[atom_id, atom_j] == dists_low_to_high[proximity]:
            neighbour = atom_j

    return neighbour


def xyz2coord(xyzs):
    """
    For a set of xyzs in the form e.g [[C, 0.0, 0.0, 0.0], ...] convert to a np array of coordinates, containing just
    just the x, y, z coordinates
    :param xyzs: List of xyzs
    :return: numpy array of coords
    """
    if isinstance(xyzs[0], list):
        return np.array([np.array(line[1:4]) for line in xyzs])
    else:
        return np.array(xyzs[1:4])


def molblock2xyzs(mol_block):
    """
    Convert an RDKit mol block to xyzs as a list of lists
    :param mol_block:
    :return: xyzs
    """

    xyzs = []

    for line in mol_block.split('\n'):
        if len(line.split()) == 16 and line.split()[0][-1].isdigit():
            x, y, z, atom_label = line.split()[:4]
            xyzs.append([atom_label, float(x), float(y), float(z)])

    return xyzs


def calc_distance_matrix(xyzs):
    """
    Calculate a distance matrix
    :param xyzs: List of xyzs
    :return:
    """

    n_atoms = len(xyzs)
    coords = xyz2coord(xyzs)
    distance_matrix = np.zeros([n_atoms, n_atoms])

    for atom_i in range(n_atoms):
        for atom_j in range(n_atoms):
            dist = np.linalg.norm(coords[atom_i] - coords[atom_j])
            distance_matrix[atom_i, atom_j] = dist

    return distance_matrix


def calc_dist(atom_i, atom_j, xyzs=None, coords=None):
    if xyzs:
        return np.linalg.norm(xyz2coord(xyzs[atom_i]) - xyz2coord(xyzs[atom_j]))
    if coords:
        return np.linalg.norm(coords[atom_i] - coords[atom_j])


def calc_midpoint(coord1, coord2):
    """
    Calculate the midpoint between two atoms
    :param coord1: Coordinate of of an atom as numpy array
    :param coord2: Coordinate of of an atom as numpy array
    :return: Coordinate of midpoint
    """
    return (coord1 + coord2) / 2.0


def is_geom_reasonable(xyzs):
    """
    For an xyz list check to ensure the geometry is sensible, before an optimisation is carried out. There should be
    no distances smaller than 0.7 Å
    :param xyzs: List of xyzs
    :return:
    """
    logger.info('Checking to see whether the geometry is reasonable')

    for n in range(len(xyzs)):
        for m in range(len(xyzs)):
            if n > m:
                line_i, line_j = xyzs[n], xyzs[m]
                dist = np.linalg.norm(np.array(line_i[1:4]) - np.array(line_j[1:4]))
                if dist < 0.8:
                    logger.warning('There is a distance < 0.8 Å. There is likely a problem with the geometry')
                    return False
                if dist > 100:
                    logger.warning('There is a distance > 100 Å. There is likely a problem with the geometry')
                    return False

    return True


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/np.linalg.norm(axis)
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def get_rot_mat_kabsch(p_matrix, q_matrix):
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


def get_centered_matrix(mat):
    """
    For a list of coordinates n.e. a n_atoms x 3 matrix as a np array translate to the center of the coordinates
    :param mat: (np.ndarray)
    :return: (np.ndarray) translated coordinates
    """
    centroid = np.average(mat, axis=0)
    return np.array([coord - centroid for coord in mat])


i = np.array([1.0, 0.0, 0.0])
j = np.array([0.0, 1.0, 0.0])
k = np.array([0.0, 0.0, 1.0])
