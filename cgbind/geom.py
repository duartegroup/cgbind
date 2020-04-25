import numpy as np
from cgbind.log import logger
from scipy.spatial import distance_matrix
from cgbind.atoms import get_atomic_mass


def calc_com(atoms):
    """
    Calculate the centre of mass for a list of xyzs

    :param atoms: (list(cgbind.atoms.Atom))
    :return: (np.ndarray) shape: 3
    """
    logger.info('Calculating centre of mass ')

    com = np.zeros(3)
    total_mass = 0.0

    for atom in atoms:
        atom_mass = get_atomic_mass(atom)
        total_mass += atom_mass

        com += atom_mass * atom.coord

    return com / total_mass


def calc_normalised_vector(coord1, coord2):

    vec = coord2 - coord1
    return vec / np.linalg.norm(vec)


def is_geom_reasonable(molecule):
    """
    For an xyz list check to ensure the geometry is sensible, before an optimisation is carried out. There should be
    no distances smaller than 0.7 Å

    :param molecule: (cgbind.molecule.BaseStruct)
    :return: (bool)
    """
    logger.info('Checking to see whether the geometry is reasonable')

    coords = molecule.get_coords()

    # Compute the distance matrix with all i,j pairs, thus add 1 to the diagonals to remove the d(ii) = 0
    # components that would otherwise result in an unreasonable geometry

    dist_mat = distance_matrix(coords, coords) + np.identity(len(coords))

    if np.min(dist_mat) < 0.8:
        logger.warning('There is a distance < 0.8 Å. There is likely a problem with the geometry')
        return False
    if np.max(dist_mat) > 1000:
        logger.warning('There is a distance > 1000 Å. There is likely a problem with the geometry')
        return False

    return True


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.

    :param axis: (np.ndarray) Unit vector in 3D to rotate around
    :param theta: (float) Angle in radians
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


def spherical_to_cart(r, theta, phi):

    return np.array([r * np.cos(theta) * np.sin(phi),
                     r * np.sin(theta) * np.sin(phi),
                     r * np.cos(phi)])


i = np.array([1.0, 0.0, 0.0])
j = np.array([0.0, 1.0, 0.0])
k = np.array([0.0, 0.0, 1.0])
