# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
from libc.math cimport sqrt, fmin
import numpy as np


cdef list get_esp_vals(int nx, int ny, int nz, double[:, :] coords, double[:] min_carts,
               double[:] charges, int n_atoms, double rx, double ry, double rz,
               double[:] dists, double[:] esp_vals,):

    cdef double min_val = 1E9
    cdef double max_val = -1E9
    cdef double min_dist
    cdef double esp_val
    cdef double vox_point[3]
    cdef double tmp
    cdef int i, j, k, n, m


    for i in range(nx):
        for j in range(ny):
            for k in range(nz):

                vox_point[0] = min_carts[0] + (<double>i + 0.5) * rx
                vox_point[1] = min_carts[1] + (<double>j + 0.5) * ry
                vox_point[2] = min_carts[2] + (<double>k + 0.5) * rz

                min_dist = 1E9
                for n in range(n_atoms):

                    # Compute the distance
                    d = 0.0
                    for m in range(3):
                        tmp = coords[n, m] - vox_point[m]
                        d += tmp * tmp

                    dists[n] = sqrt(d)

                    if dists[n] < min_dist:
                        min_dist = dists[n]

                # ESP_j = q_i / r_ij
                esp_val = 0.0
                for n in range(n_atoms):
                    esp_val += charges[n] / dists[n]

                # Min and max values are on *roughly* the van der Walls surface i.e. a distance r from the closest atom
                # where r is the VdW radius of the closest atom. Taken as 1.5 -> 2.0 Å here
                if 1.5 < min_dist < 2.0:
                    if esp_val < min_val:
                        min_val = esp_val
                if 1.5 < min_dist < 2.0:
                    if esp_val > min_val:
                        max_val = esp_val

                esp_vals[i*nx*ny + j * ny + k] = esp_val

    cdef list results = [esp_vals, min_val, max_val]
    return results



def get_cube_lines(py_nx, py_ny, py_nz, py_coords, py_min_carts, py_charges, py_vox_size):
    """
    Generate the value part of the electrostatic potential .cube file

    :param py_nx: (int) Number of points in x
    :param py_ny: (int) Number of points in y
    :param py_nz: (int) Number of points in z
    :param py_coords: (list(np.ndarray)) Coordinates of the atoms in bohr (Nx3)
    :param py_min_carts: (np.ndarray) Bottom left corner of the strucure to start
                         generating the ESP from (length 3)
    :param py_charges: (np.ndarray) List of partial atomic charges (length N)
    """

    cdef int nx = py_nx
    cdef int ny = py_ny
    cdef int nz = py_nz

    cdef double rx = py_vox_size[0] / py_nx
    cdef double ry = py_vox_size[1] / py_ny
    cdef double rz = py_vox_size[2] / py_nz

    cdef double[:, :] coords = np.array(py_coords, dtype='f8')
    cdef double[:] min_carts = np.array(py_min_carts, dtype='f8')
    cdef double[:] charges = np.array(py_charges, dtype='f8')
    cdef int n_atoms = len(py_coords)

    cdef double[:] dists = np.zeros(len(coords), dtype='f8')
    cdef double[:] esp_vals = np.zeros(py_nx*py_ny*py_nz, dtype='f8')

    res = get_esp_vals(nx, ny, nz, coords, min_carts, charges, n_atoms, rx, ry,
                       rz, dists, esp_vals)
    res[0] = esp_vals

    # Construct the ESP file lines from the flat list of values. Fast so no need for c implementation
    line_list = []
    for i in range(py_nx):
        for j in range(py_ny):
            for k in range(py_nz):
                esp_val = esp_vals[i*py_nx*py_ny + j * py_ny + k]
                line_list.append(f'{esp_val:>15.5E}  ')
                if k % 6 == 5:
                    line_list.append('\n')

            line_list.append('\n')


    # Cube file lines, min ESP value and max ESP value close to the VdW radius
    return line_list, res[1], res[2]
