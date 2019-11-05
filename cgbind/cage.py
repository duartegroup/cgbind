from copy import deepcopy
import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import basinhopping, minimize
from cgbind.x_motifs import get_shifted_template_x_motif_coords
from cgbind.build import get_fitted_linker_coords
from cgbind.log import logger
from cgbind.input_output import print_output
from cgbind.atoms import get_vdw_radii
from cgbind.geom import is_geom_reasonable
from cgbind.geom import xyz2coord
from cgbind.geom import spherical_to_cart
from cgbind import calculations
from cgbind.input_output import xyzs2xyzfile


def get_max_sphere_negative_radius(theta_and_phi, r, cage_coords):
    """
    Get the maximum sphere radius that is possible at a point defined by the spherical polar coordinates theta, phi and
    r. This amounts to finding the minimum pairwise distance between the point and the rest of the cage. The negative
    radius is returned as it will be fed into scipy.optmize.minimise

    :param theta_and_phi: (list(float))
    :param r: (float)
    :param cage_coords: (np.ndarray) n_atoms x 3
    :return: (float)
    """

    theta, phi = theta_and_phi
    # Convert the point in spherical polars to Cartesian so the distances to the rest of the cage can be calculated
    # nneds to be a 1 x 3 matrix to use cdist
    point = np.array([spherical_to_cart(r=r, theta=theta, phi=phi)])

    return -np.min(cdist(point, cage_coords))


class Cage(object):

    def _is_linker_reasonable(self, linker):

        if linker is None:
            logger.error(f'Linker was None. Cannot build {self.name}')
            return False

        if linker.xyzs is None or linker.arch is None or linker.name is None:
            logger.error(f'Linker doesn\'t have all the required attributes. Cannot build {self.name}')
            return False

        return True

    def get_centroid(self):
        metal_coords = [xyz2coord(xyz) for xyz in self.xyzs if self.metal in xyz]
        return np.average(metal_coords, axis=0)

    def print_xyzfile(self, force=False):
        if self.reasonable_geometry or force:
            xyzs2xyzfile(xyzs=self.xyzs, basename=self.name)

    def get_metal_atom_ids(self):
        logger.info(f'Getting metal_label atom ids with label {self.metal}')
        try:
            return [i for i in range(len(self.xyzs)) if self.xyzs[i][0] == self.metal]
        except TypeError or IndexError or AttributeError:
            logger.error('Could not get metal_label atom ids. Returning None')
            return None

    def get_cavity_vol(self):
        """
        For a cage extract the cavity volume defined as the volume of the largest sphere, centered on the M-M midpoint
        that may be constructed while r < r(midpoint--closest atom)
        :return: Cavity volume in Å^3
        """
        logger.info('Calculating maximum enclosed sphere')

        min_centriod_atom_dist = 999.9
        centroid, min_atom_dist_id = None, None

        try:
            centroid = self.get_centroid()
            if centroid is None:
                logger.error('Could not find the cage centroid. Returning 0.0')
                return 0.0

            # Compute the smallest distance to the centroid
            for i in range(len(self.xyzs)):
                dist = np.linalg.norm(xyz2coord(self.xyzs[i]) - centroid)
                if dist < min_centriod_atom_dist:
                    min_centriod_atom_dist = dist
                    min_atom_dist_id = i

        except TypeError or ValueError or AttributeError:
            pass

        if min_atom_dist_id is not None:
            vdv_radii = get_vdw_radii(atom_label=self.xyzs[min_atom_dist_id][0])
            # V = 4/3 π r^3, where r is ithe centroid -> closest atom distance, minus it's VdW volume
            return (4.0 / 3.0) * np.pi * (min_centriod_atom_dist - vdv_radii)**3

        else:
            logger.error('Could not calculate the cavity volume. Returning 0.0')
            return 0.0

    def get_m_m_dist(self):
        """
        For a cage calculate the M-M distance
        :return: Distance in Å
        """
        try:
            m_m_dists = []
            for m_id_i in range(len(self.m_ids)):
                for m_id_j in range(len(self.m_ids)):
                    if m_id_i > m_id_j:
                        dist = np.linalg.norm(xyz2coord(self.xyzs[self.m_ids[m_id_i]]) -
                                              xyz2coord(self.xyzs[self.m_ids[m_id_j]]))
                        m_m_dists.append(dist)

            if len(m_m_dists) > 0:
                return np.average(np.array(m_m_dists))
            else:
                logger.error('Could not find any metal_label atoms')

        except TypeError or ValueError or AttributeError:
            logger.error('Could not calculate the M-M distance. Returning 0.0')

        return 0.0

    def get_n_rot_bonds(self):
        try:
            return sum([linker.n_rot_bonds for linker in self.linkers])
        except TypeError:
            return None

    def get_h_bond_donors(self):
        try:
            return sum([linker.n_h_donors for linker in self.linkers])
        except TypeError:
            return None

    def get_max_escape_sphere(self, basinh=False, max_dist_from_metals=10):
        """
        Get the maximum radius of a sphere that can escape from the centroid of the cage – will iterate through all
        theta/phi

        :param basinh: (bool) Find the true maximum escape sphere by basin hopping on the surface
        :param max_dist_from_metals: (float) Distance in Å on top of the average M-M distance that will be used for the
                                             search for the maximum escape sphere
        :return:
        """
        logger.info('Getting the volume of the largest sphere that can escape from the cavity')

        max_sphere_escape_r = 99999999999999.9
        avg_m_m_dist = self.get_m_m_dist()
        cage_coords = xyz2coord(self.xyzs) - self.get_centroid()

        # For a distance from the origin (the cage centroid) calculate the largest sphere possible without hitting atoms
        opt_theta_phi, opt_r = np.zeros(2), 0.0
        for r in np.linspace(0.0, avg_m_m_dist + max_dist_from_metals, 20):
            if basinh:
                opt = basinhopping(get_max_sphere_negative_radius, x0=opt_theta_phi, stepsize=1.0, niter=5,
                                   minimizer_kwargs={'args': (r, cage_coords), 'method': 'BFGS'})
            else:
                opt = minimize(get_max_sphere_negative_radius, x0=opt_theta_phi, args= (r, cage_coords), method='BFGS')

            opt_theta_phi = opt.x

            # This is the correct way round because we want the largest sphere that CAN escape
            if -opt.fun < max_sphere_escape_r:
                max_sphere_escape_r = -opt.fun
                opt_r = r

        print(max_sphere_escape_r)

        # Get the atom id that the max escape sphere hits into
        sphere_point = spherical_to_cart(r=opt_r, theta=opt_theta_phi[0], phi=opt_theta_phi[1])
        atom_id = np.argmin([np.linalg.norm(coord - sphere_point) for coord in cage_coords])

        radius = max_sphere_escape_r - get_vdw_radii(atom_label=self.xyzs[atom_id][0])
        logger.info(f'Radius of largest sphere that can escape from the cavity = {radius}')

        return (4.0 / 3.0) * np.pi * radius**3

    def singlepoint(self, method, keywords, n_cores=1, max_core_mb=1000):
        return calculations.singlepoint(self, method, keywords, n_cores, max_core_mb)

    def optimise(self, method, keywords, n_cores=1, max_core_mb=1000):
        return calculations.optimise(self, method, keywords, n_cores, max_core_mb)

    def _calc_charge(self):
        logger.info('Calculating the charge on the metallocage')
        self.charge = self.arch.n_metals * self.metal_charge + np.sum(np.array([linker.charge for linker in self.linkers]))
        return None

    def _init_homoleptic_cage(self, linker):
        logger.info(f'Initialising a homoleptic cage')

        if not self._is_linker_reasonable(linker):
            logger.error('Linker was not reasonable')
            return

        self.name = 'cage_' + linker.name
        self.arch = linker.arch
        self.dr = linker.dr
        self.linkers = [linker for _ in range(linker.arch.n_linkers)]
        self.cage_template = linker.cage_template

        return

    def _init_heteroleptic_cage(self, linkers):
        logger.info(f'Initialising a heteroleptic cage')

        if not all([self._is_linker_reasonable(linker) for linker in linkers]):
            logger.error('Not all linkers were reasonable')
            return

        if not all([linker.arch.name == linkers[0].arch.name for linker in linkers]):
            logger.error('Linkers had different architectures, not building a cage')
            return

        self.name = 'cage_' + '_'.join([linker.name for linker in linkers])
        self.arch = linkers[0].arch
        self.linkers = linkers
        self.cage_template = linkers[0].cage_template

        logger.warning('Hetroleptic cages will have the average dr of all linkers')
        self.dr = np.average(np.array([linker.dr for linker in linkers]))

        return

    def build(self):
        logger.info('Building a cage geometry')

        if self.dr is None:
            logger.error('Cannot build a cage dr was None')
            return None

        xyzs = []
        # Add the metals from the template shifted by dr
        for metal in self.cage_template.metals:
            metal_coord = self.dr * metal.shift_vec / np.linalg.norm(metal.shift_vec) + metal.coord
            xyzs.append([self.metal] + metal_coord.tolist())

        # Add the linkers by shifting the x_motifs in each linker templates by dr and finding the best rot matrix
        for i, template_linker in enumerate(self.cage_template.linkers):

            new_linker = deepcopy(self.linkers[i])
            shifted_coords = get_shifted_template_x_motif_coords(linker_template=template_linker, dr=self.dr)
            x_coords = [new_linker.coords[atom_id] for motif in new_linker.x_motifs for atom_id in motif.atom_ids]

            linker_coords  = get_fitted_linker_coords(linker=self.linkers[i], template_x_coords=shifted_coords,
                                                      coords_to_fit=x_coords, current_xyzs=xyzs)

            xyzs += [[new_linker.xyzs[i][0]] + linker_coords[i].tolist() for i in range(new_linker.n_atoms)]

        if len(xyzs) != self.arch.n_metals + np.sum(np.array([linker.n_atoms for linker in self.linkers])):
            logger.error('Failed to build a cage')
            return None

        return xyzs

    def __init__(self, linker=None, metal=None, metal_charge=0, linkers=None, solvent=None, mult=1):
        """
        Initialise a cage object
        :param linker: Linker object
        :param metal: (str)
        :param metal_charge: (int)
        """
        logger.info(f'Initialising a Cage object')

        self.metal = metal
        self.mult = mult
        self.solvent = solvent
        self.name = 'cage'                                                           # Will be overwritten in _init_cage
        self.linkers = None
        self.dr = None
        self.arch = None
        self.cage_template = None

        self.metal_charge = int(metal_charge)
        self.charge = None
        self.energy, self.xyzs, self.m_ids = None, None, None

        self.reasonable_geometry = False

        if linker is not None:
            self._init_homoleptic_cage(linker)

        elif linkers is not None:
            self._init_heteroleptic_cage(linkers)

        else:
            logger.error('Could not generate a cage object without either a linker or set of linkers')
            return

        if self.linkers is None:
            logger.error('Cannot build a cage with linkers as None')
            return

        self._calc_charge()

        self.reasonable_geometry = True
        self.xyzs = self.build()

        if self.xyzs is None:
            self.reasonable_geometry = False
            print_output('Cage build for', self.name, 'Failed')
            return

        self.m_ids = self.get_metal_atom_ids()
        self.reasonable_geometry = is_geom_reasonable(self.xyzs)

        print_output('Cage', self.name, 'Built')
