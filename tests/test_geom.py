from cgbind import geom
import numpy as np

xyz_list = [['H', 0.0, 0.0, 0.0], ['H', 1.0, 0.0, 0.0]]


def test_xyz2coord():

    coord_list = geom.xyz2coord(xyz_list)

    assert type(coord_list) == np.ndarray
    assert type(coord_list[0]) == np.ndarray
    assert 0.99 < np.linalg.norm(coord_list[0] - coord_list[1]) < 1.001

    xyz_line = ['H', 0.0, 0.0, 0.0]
    coord = geom.xyz2coord(xyz_line)

    assert type(coord) == np.ndarray
    assert len(coord) == 3


def test_distance_matrix():

    distance_matix = geom.calc_distance_matrix(xyz_list)
    assert distance_matix.shape == (2, 2)
    assert distance_matix[0, 0] == 0.0


def test_rot_matix():

    axis = np.array([0.0, 0.0, 1.0])
    theta = np.pi                       # angle in radians

    rot_matix = geom.rotation_matrix(axis, theta)
    point = np.array([1.0, 1.0, 1.0])
    rot_point = np.matmul(rot_matix, point)

    assert rot_matix.shape == (3, 3)
    assert -1.001 < rot_matix[0, 0] < -0.999
    assert -1.001 < rot_matix[1, 1] < -0.999
    assert 0.999 < rot_matix[2, 2] < 1.001
    assert -0.001 < rot_matix[0, 1] < 0.001

    assert -1.001 < rot_point[0] < -0.999
    assert -1.001 < rot_point[1] < -0.999
    assert 0.999 < rot_point[2] < 1.001


def test_com():
    com = geom.calc_com(xyzs=xyz_list)
    ideal_com = np.array([0.5, 0.0, 0.0])
    assert np.abs(np.average(com - ideal_com)) < 1E-5


def test_normed_vector():

    coord1 = np.array([0.0, 0.0, 0.0])
    coord2 = np.array([2.0, 0.0, 0.0])

    ideal_normed_vector = np.array([1.0, 0.0, 0.0])
    normed_vector = geom.calc_normalised_vector(coord1, coord2)

    assert np.abs(np.average(ideal_normed_vector - normed_vector)) < 1E-5


def test_closest_atom_id():
    bonded_h2 = [['H', 0.0, 0.0, 0.0], ['H', 0.8, 0.0, 0.0]]
    assert geom.get_closest_bonded_atom_id(xyzs=bonded_h2, atom_id=0) == 1


def test_dist_calc():
    ideal_length = 1.0
    assert np.abs(geom.calc_dist(xyzs=xyz_list, atom_i=0, atom_j=1) - ideal_length) < 1E-5


def test_reasonable_geom():
    xyzs = xyz_list.copy()
    xyzs.append(['H', 0.5, 0.0, 0.0])
    assert geom.is_geom_reasonable(xyzs=xyzs) is False
    assert geom.is_geom_reasonable(xyzs=xyz_list) is True


def test_neighbour():
    xyzs = xyz_list.copy()
    xyzs.append(['H', 0.5, 0.0, 0.0])
    dist_matrix = geom.calc_distance_matrix(xyzs)
    assert geom.get_neighbour(atom_id=0, distance_matrix=dist_matrix, proximity=1) == 2
    assert geom.get_neighbour(atom_id=0, distance_matrix=dist_matrix, proximity=2) == 1
