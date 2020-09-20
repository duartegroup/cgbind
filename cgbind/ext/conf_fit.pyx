# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
from libc.math cimport cos, sin, sqrt
import numpy as np

cdef double [:, :] update_rotation_matrix(double [:, :] rot_mat,
                                          double [:] axis,
                                          double theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.

    :param rot_mat:
    :param axis: (np.ndarray) Unit vector in 3D to rotate around
    :param theta: (float) Angle in radians
    """
    cdef double norm = sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
    cdef double sin_theta = sin(theta/2.0)

    cdef double a = cos(theta/2.0)
    cdef double b = -axis[0] * sin_theta / norm
    cdef double c = -axis[1] * sin_theta / norm
    cdef double d = -axis[2] * sin_theta / norm

    rot_mat[0, 0] = a*a+b*b-c*c-d*d
    rot_mat[0, 1] = 2.0*(b*c+a*d)
    rot_mat[0, 2] = 2.0*(b*d-a*c)

    rot_mat[1, 0] = 2.0*(b*c-a*d)
    rot_mat[1, 1] = a*a+c*c-b*b-d*d
    rot_mat[1, 2] = 2.0*(c*d+a*b)

    rot_mat[2, 0] = 2.0*(b*d+a*c)
    rot_mat[2, 1] = 2.0*(c*d-a*b)
    rot_mat[2, 2] = a*a+d*d-b*b-c*c

    return rot_mat


cdef void coords_rotated_dihedral(double [:, :] coords,
                                  double[:] axis,
                                  double[:] origin,
                                  double theta,
                                  int[:] rot_idxs,
                                  double [:, :] rot_mat,
                                  int n_atoms):
    """
    Generate coordinates by rotation around a bond

    :param coords:
    :param axis:
    :param origin:
    :param theta: (float)
    :param rot_idxs:
    :param rot_mat:
    :param n_atoms:
    :return:
    """
    rot_mat = update_rotation_matrix(rot_mat,
                                     axis=axis,
                                     theta=theta)
    cdef int i, j, k
    cdef double x, y, z

    # Rotate the portion of the structure to rotate about the bonding axis
    for i in range(n_atoms):

        if rot_idxs[i] == 1:
            for k in range(3):
                coords[i][k] -= origin[k]

            # Apply the rotation
            x = coords[i][0]
            y = coords[i][1]
            z = coords[i][2]

            for k in range(3):
                coords[i][k] = (rot_mat[k, 0] * x +
                                rot_mat[k, 1] * y +
                                rot_mat[k, 2] * z)

            for k in range(3):
                coords[i][k] += origin[k]
    return


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


def rotate_coords(py_coords, py_axis, py_theta, py_origin, py_rot_idxs):
    """
    Rotate some atoms by an amount theta radians in an axis

    :param py_axis:
    :param py_theta:
    :param py_origin:
    :param py_coords:
    :param py_rot_idxs:
    :return:
    """
    cdef double [:] axis = np.array(py_axis, dtype='f8')
    cdef double [:] origin = np.array(py_origin, dtype='f8')
    cdef double theta = py_theta
    cdef int [:] rot_idxs = np.zeros(len(py_coords), dtype=np.int32)

    cdef double [:, :] coords = np.array(py_coords, dtype='f8')
    cdef double [:, :] rot_mat = np.zeros(shape=(3, 3), dtype='f8')

    cdef int n_atoms = len(py_coords)

    for i in range(len(py_coords)):
        if i in py_rot_idxs:
            rot_idxs[i] = 1

    coords_rotated_dihedral(coords=coords,
                            axis=axis,
                            origin=origin,
                            theta=theta,
                            rot_idxs=rot_idxs,
                            rot_mat=rot_mat,
                            n_atoms=n_atoms)

    return np.asarray(coords)


def opt_coords(thetas_list, bonds_and_rot_idxs, py_coords, py_template_coords,
               fit_idxs):
    """
    Calculate the list of rotations about dihedral bonds that minimise the
    fitting cost betwee
    n the coordinates of the X-motif atoms and the template
    coordinates. The list of possibilities will be enumerated exhaustively.
    Return the optimal fitting coordinates

    :param thetas_list: (list(list(float)))
    :param bonds_and_rot_idxs: (dict) Keys of bonds (tuple(int)) and values of
                              the atom indexes right of the bond to rotate
    :param py_coords: (np.ndarray)
    :param py_template_coords: (np.ndarray)
    :param fit_idxs: (list(int))
    :return: (list(float))
    """
    # A set of thetas should be the same as the number of bonds that can be
    # rotated
    assert len(thetas_list[0]) == len(bonds_and_rot_idxs)

    cdef double [:, :] thetas = np.array(thetas_list, dtype='f8')
    cdef int n_thetas = len(thetas_list)

    # Matrix of bond indexes that will be rotated around
    cdef int [:, :] bonds = np.zeros(shape=(len(bonds_and_rot_idxs), 2),
                                     dtype=np.int32)
    cdef int n_bonds = len(bonds_and_rot_idxs)

    for n, bond in enumerate(bonds_and_rot_idxs):
        for m in range(2):
            bonds[n, m] = bond[m]

    # Matrix of indexes that can be rotated set with 1 if it is rotated
    # and zero otherwise
    cdef int [:, :] rot_idxs = np.zeros(shape=(len(bonds_and_rot_idxs),
                                               len(py_coords)), dtype=np.int32)
    cdef int n_atoms = len(py_coords)

    for n, idxs in enumerate(bonds_and_rot_idxs):
        for m in range(len(py_coords)):
            if m in idxs:
                rot_idxs[n, m] = 1

    cdef double [:, :] coords = np.array(py_coords, dtype='f8')
    cdef double [:, :] rot_coords = np.array(py_coords, dtype='f8')

    cdef double [:, :] t_coords = np.array(py_template_coords, dtype='f8')

    cdef double [:, :] rot_mat = np.zeros(shape=(3, 3), dtype='f8')
    cdef double [:] axis = np.zeros(3, dtype='f8')
    cdef double [:] origin = np.zeros(3, dtype='f8')

    cdef int l_idx, r_idx
    cdef int i, j, k, l

    min_cost, opt_coords = None, None

    for i in range(n_thetas):

        # Reset the coordinates
        for k in range(n_atoms):
            for l in range(3):
                rot_coords[k, l] = coords[k, l]

        for j in range(n_bonds):

            l_idx = bonds[j][0]
            r_idx = bonds[j][1]

            # Set the axis and origin
            for k in range(3):
                axis[k] = coords[l_idx][k] - coords[r_idx][k]
                origin[k] = coords[l_idx][k]

            # And rotate around this dihedral
            coords_rotated_dihedral(coords=rot_coords,
                                    axis=axis,
                                    origin=coords[l_idx],
                                    theta=thetas[i, j],
                                    rot_idxs=rot_idxs[j],
                                    rot_mat=rot_mat,
                                    n_atoms=n_atoms)

        rot_coords = np.asarray(rot_coords)
        cost = fit_cost(np.array([rot_coords[n, :] for n in range(len(py_coords)) if n in fit_idxs]), t_coords)
        print(cost)

        if min_cost is None or cost < min_cost:
            min_cost = cost
            py_coords = np.array([rot_coords[n, :] for n in range(len(py_coords))])

    return py_coords
