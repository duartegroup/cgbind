from cgbind import geom
import numpy as np
from cgbind.atoms import Atom
from cgbind.molecule import Molecule

h2_atoms = [Atom('H', 0.0, 0.0, 0.0), Atom('H', 1.0, 0.0, 0.0)]


def test_com():
    com = geom.calc_com(atoms=h2_atoms)
    ideal_com = np.array([0.5, 0.0, 0.0])
    assert np.abs(np.average(com - ideal_com)) < 1E-5


def test_normed_vector():

    coord1 = np.array([0.0, 0.0, 0.0])
    coord2 = np.array([2.0, 0.0, 0.0])

    ideal_normed_vector = np.array([1.0, 0.0, 0.0])
    normed_vector = geom.calc_normalised_vector(coord1, coord2)

    assert np.abs(np.average(ideal_normed_vector - normed_vector)) < 1E-5


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


def test_reasonable_geom():

    tmp = Molecule()
    tmp.set_atoms([Atom('H', 0.0, 0.0, 0.0), Atom('H', 0.1, 0.0, 0.0)])
    assert geom.is_geom_reasonable(tmp) is False
    tmp.set_atoms([Atom('H', 0.0, 0.0, 0.0), Atom('H', 1.0, 0.0, 0.0)])
    assert geom.is_geom_reasonable(tmp) is True
    tmp.set_atoms([Atom('H', 0.0, 0.0, 0.0), Atom('H', 1001, 0.0, 0.0)])
    assert geom.is_geom_reasonable(tmp) is False


def test_rot_kabsch():

    p_mat = np.array([[ 5.84606667, -0.7098, -0.42831667],
                      [ 6.94516667, -0.0456, -0.80221667],
                      [ 7.36736667, 1.0454, -0.16211667],
                      [-5.88033333, -0.835, 0.01198333],
                      [-7.05143333, -0.1629, 0.18438333],
                      [-7.22683333, 0.7079, 1.19628333]])

    q_mat = np.array([[-2.19739193, -5.39769871, -0.81494868],
                     [-3.17969193, -6.24319871, -1.15854868],
                     [-2.92989193, -7.56459871, -1.25704868],
                     [2.50739193, 5.25459871, 0.97654868],
                     [2.36369193, 6.58479871, 1.00314868],
                     [3.43589193, 7.36609871, 1.25084868]])

    rot_mat = np.array([[-0.34283885, 0.49631223, 0.79758115],
                        [-0.933562, -0.08554556, -0.34805739],
                        [-0.10451561, -0.86391905, 0.49266658]])

    assert np.sum(geom.get_rot_mat_kabsch(p_matrix=p_mat, q_matrix=q_mat) - rot_mat) < 0.01


def test_get_centered_matrix():

    mat = np.array([[0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0]])

    assert -0.51 < geom.get_centered_matrix(mat)[0][0] < -0.49


def test_spherical_to_cart():

    assert np.sum(geom.spherical_to_cart(r=1.0, theta=0.0, phi=0.0) - np.array([0.0, 0.0, 1.0])) < 0.01
