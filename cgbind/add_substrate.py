from cgbind.log import logger
from copy import deepcopy
import numpy as np
from cgbind.constants import Constants
from rdkit.Chem import AllChem
from scipy.optimize import minimize, Bounds
from scipy.spatial import distance_matrix
from cgbind import geom
from cgbind.atoms import get_vdw_radii
from cgbind.geom import rotation_matrix
from cgbind.geom import calc_com
from cgbind.utils import copy_func


def cage_subst_repulsion_func(cage, substrate, cage_coords, subst_coords, with_attraction=True):
    """
    Determine the energy using two-body atom-atom repulsion derived from noble gas dimers where

    V_rep(r) = exp(- r/b + a)

    where a and b are parameters determined by the atom pairs. Parameters are suitable to generate V_rep in kcal mol-1

    :param cage: (Cage object)
    :param substrate: (Substrate object)
    :param cage_coords: (list(np.ndarray)) Cage coordinates
    :param subst_coords: (list(np.ndarray)) Substrate coordinates
    :param with_attraction: (bool) do or don't return the energy with a constant attractive term based on the number of
                                   substrate atoms in the structure

    :return: energy: (float) Potential energy (V_rep) in kcal mol-1
    """

    dist_mat = distance_matrix(cage_coords, subst_coords)

    # Matrix with the pairwise additions of the vdW radii
    sum_vdw_radii = np.add.outer(np.array(cage.vdw_radii),
                                 np.array(substrate.vdw_radii))

    # Magic numbers derived from fitting potentials to noble gas dimers and plotting against the sum of vdw radii
    b_mat = 0.083214 * sum_vdw_radii - 0.003768
    a_mat = 11.576415 * (0.175541 * sum_vdw_radii + 0.316642)
    exponent_mat = -(dist_mat / b_mat) + a_mat

    energy_mat = np.exp(exponent_mat)
    energy = np.sum(energy_mat)

    #      E is negative for favourable binding but this is a purely repulsive function so subtract a number..
    #      which is determined from the best classifier for 102 binding affinities (see cgbind paper) 0.4 kcal mol-1
    if with_attraction:
        return energy - 0.4 * substrate.n_atoms

    return energy


def cage_subst_repulsion_and_electrostatic_func(cage, substrate, cage_coords, subst_coords):
    """
    Determine the energy of adding a substrate to a cage based on V_rep + V_att where the attractive term
    is electrostatic and uses the sum of q_i q_j / r_ij interaction energies where q_i is the partial atomic
    charge on atom i.

    :param cage: (Cage object)
    :param substrate: (Substrate object)
    :param cage_coords: (list(np.ndarray)) Cage coordinates
    :param subst_coords: (list(np.ndarray)) Substrate coordinates
    :return:
    """
    # Calculate the distance matrix in Bohr (a0) so the energies are in au
    dist_mat = Constants.ang2a0 * distance_matrix(cage_coords, subst_coords)

    # Charges are already in units of e
    prod_charge_mat = np.outer(cage.charges, substrate.charges)

    # Compute the pairwise iteration energies as V = q1 q2 / r in atomic units
    energy_mat = prod_charge_mat / dist_mat
    electrostatic_energy = Constants.ha2kcalmol * np.sum(energy_mat)

    repulsive_energy = cage_subst_repulsion_func(cage, substrate, cage_coords, subst_coords)

    return electrostatic_energy + repulsive_energy


def add_substrate_com(cagesubt):
    """
    Add a substrate the centre of a cage defined by its centre of mass (com) will minimise the energy with respect to
    rotation of the substrate and the substrate conformer using cagesubt.energy_func. Will rotate cagesubt.n_init_geom
    times and use cagesubt.n_subst_confs number of substrate conformers

    :param cagesubt: (CageSubstrateComplex object)

    :return: xyzs: (list(list))
    """
    logger.info(f'Adding substrate to the cage COM and minimising the energy with {cagesubt.energy_func.__name__}')

    # Minimum energy initialisation and the x parameter array (angles to rotate about the x, y, z axes)
    min_energy, curr_x = 9999999999.9, np.zeros(3)

    # Optimum (minimum energy) conformer
    best_coords = None

    c, s = cagesubt.cage, cagesubt.substrate
    cage_coords = get_centered_cage_coords(c)
    c.vdw_radii = [get_vdw_radii(atom) for atom in c.atoms]

    if cagesubt.n_subst_confs > 1:
        try:
            s.gen_confs(n_confs=cagesubt.n_subst_confs)
        except (ValueError, RuntimeError):
            logger.error('Could not generate substrate conformers')
            return None

    for i, substrate in enumerate(s.conformers):
        subst_coords = get_centered_substrate_coords(substrate)
        s.vdw_radii = [get_vdw_radii(atom) for atom in s.atoms]
        if s.mol_obj is not None:
            s.volume = AllChem.ComputeMolVolume(s.mol_obj, confId=i)

        for _ in range(cagesubt.n_init_geom):
            rot_angles = 2.0 * np.pi * np.random.rand(3)        # rand generates in [0, 1] so multiply with

            # Minimise the energy with a BFGS minimiser supporting bounds on the values (rotation is periodic)
            result = minimize(get_energy, x0=np.array(rot_angles),
                              args=(c, s, cagesubt.energy_func, cage_coords, subst_coords), method='L-BFGS-B',
                              bounds=Bounds(lb=0.0, ub=2*np.pi), tol=0.01)

            energy = result.fun
            logger.info(f'Energy = {energy:.4f}')

            if energy < min_energy:
                min_energy = energy
                best_coords = get_rotated_subst_coords(result.x, subst_coords)

    logger.info(f'Min energy = {min_energy:.4f} kcal mol-1')
    cagesubt.binding_energy_kcal = min_energy

    if best_coords is not None:
        s.set_atoms(coords=best_coords)
        c.set_atoms(coords=cage_coords)

        return c.atoms + s.atoms

    else:
        return None


def get_centered_cage_coords(cage):
    """Get the cage coordinates that had been translated to the cage centroid"""

    cage_coords = cage.get_coords()
    centroid = cage.get_centroid()

    return np.array([coord - centroid for coord in cage_coords])


def get_centered_substrate_coords(substrate):
    """Get the substrate coordinates that have been translated to its center of mass"""

    substrate.centre()
    return substrate.get_coords()


def cat_cage_subst_coords(cage, substrate, cage_coords, substrate_coords):
    """
    Concatenate some coordinates into a set of xyzs by adding back the atom labels from the original xyzs

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
    """Get substrate coordinates that have been rotated by x[0] radians in the x axis etc."""

    x_rot, y_rot, z_rot = x

    rot_matrix = np.identity(3)
    rot_matrix = np.matmul(rot_matrix, rotation_matrix(axis=geom.i, theta=x_rot))
    rot_matrix = np.matmul(rot_matrix, rotation_matrix(axis=geom.j, theta=y_rot))
    rot_matrix = np.matmul(rot_matrix, rotation_matrix(axis=geom.k, theta=z_rot))

    return np.array([np.matmul(rot_matrix, coord) for coord in deepcopy(subst_coords)])


def get_energy(x, cage, substrate, energy_func, cage_coords, subst_coords):
    """
    Calculate the energy in kcal mol-1 for a particular x, which contains the rotations in x, y, z cartesian directions
    """
    rot_substrate_coords = get_rotated_subst_coords(x, subst_coords)
    energy = energy_func(cage, substrate, cage_coords, rot_substrate_coords)

    return energy


cage_subst_repulsion_func.__name__ = 'repulsion'
cage_subst_repulsion_and_electrostatic_func.__name__ = 'electrostatic'
cage_subst_repulsion_and_electrostatic_func_est = copy_func(cage_subst_repulsion_and_electrostatic_func)
cage_subst_repulsion_and_electrostatic_func_est.__name__ = 'electrostatic_fast'


energy_funcs = [cage_subst_repulsion_func,
                cage_subst_repulsion_and_electrostatic_func,
                cage_subst_repulsion_and_electrostatic_func_est]
