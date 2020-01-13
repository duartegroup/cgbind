import numpy as np
from copy import deepcopy
from scipy.spatial.distance import cdist
from scipy.optimize import basinhopping, minimize
from cgbind.molecule import BaseStruct
from cgbind.calculations import get_charges
from cgbind.build import get_linker_xyzs_to_add_and_cost
from cgbind.log import logger
from cgbind.input_output import print_output
from cgbind.atoms import get_vdw_radii
from cgbind.geom import is_geom_reasonable
from cgbind.geom import xyz2coord
from cgbind.geom import spherical_to_cart
from cgbind.esp import get_esp_cube_lines


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


class Cage(BaseStruct):

    def get_centroid(self):
        """
        Get the centroid of a metallocage. Defined as the midpoint between all self.metal atoms in the structure

        :return: (np.ndarray) Centroid coordinate (x, y, z)
        """

        metal_coords = [xyz2coord(xyz) for xyz in self.xyzs if self.metal in xyz]
        return np.average(metal_coords, axis=0)

    def get_esp_cube(self):
        """
        Get the electrostatic potential (ESP) in a Gaussian .cube format by calculating partial atomic charges using
        XTB (tested with v. 6.2). Calls self.get_charges() and depends on self.xyzs

        :return: (list) .cube file lines
        """
        esp_lines = get_esp_cube_lines(charges=self.get_charges(), xyzs=self.xyzs)
        return esp_lines

    def print_esp_cube_file(self):
        """
        Print an electrostatic potential (ESP) .cube file. Prints the lines from self.get_esp_cube()

        :return: None
        """

        cube_file_lines = self.get_esp_cube()

        if len(cube_file_lines) == 0:
            logger.error('Could not generate cube')
            return None

        with open(self.name + '_esp.cube', 'w') as cube_file:
            [print(line, end='', file=cube_file) for line in cube_file_lines]

        return None

    def get_charges(self, estimate=False, guess=False):
        """
        Get the partial atomic charges on the cage either using XTB or estimate using no polarisation i.e. the metals
        retain their full charge and the linker charges are estimated using the Gasteiger scheme in RDKit

        :param estimate: (bool)
        :param guess: (bool) Guess the charges based on the electronegativity
        :return: (function) calculations.get_charges(self)
        """
        if estimate or guess:

            charges = [self.metal_charge] * self.arch.n_metals
            for linker in self.linkers:
                linker_charges = linker.get_charges(estimate=estimate, guess=guess)
                charges += linker_charges

            return charges

        return get_charges(self)

    def get_metal_atom_ids(self):
        """
        Get the atom ids of the self.metal atoms in the xyzs

        :return: (list(int))
        """

        logger.info(f'Getting metal_label atom ids with label {self.metal}')
        try:
            return [i for i in range(len(self.xyzs)) if self.xyzs[i][0] == self.metal]
        except TypeError or IndexError or AttributeError:
            logger.error('Could not get metal_label atom ids. Returning None')
            return None

    def get_cavity_vol(self):
        """
        For a cage extract the cavity volume defined as the volume of the largest sphere, centered on the cage centroid
        that may be constructed while r < r(midpoint--closest atom)

        :return: (float) Cavity volume in Å^3
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
            # V = 4/3 π r^3, where r is the centroid -> closest atom distance, minus it's VdW volume
            return (4.0 / 3.0) * np.pi * (min_centriod_atom_dist - vdv_radii)**3

        else:
            logger.error('Could not calculate the cavity volume. Returning 0.0')
            return 0.0

    def get_m_m_dist(self):
        """
        For a cage calculate the average M-M distance

        :return: (float) Distance in Å
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

    def get_num_rot_bonds(self):
        """
        Get the number of rotatable bonds in a metallocage

        :return: (int)
        """
        try:
            return sum([linker.n_rot_bonds for linker in self.linkers])
        except TypeError:
            return None

    def get_num_h_bond_donors(self):
        """
        Get the number of hydrogen bond donors in a metallocage

        :return: (int)
        """
        try:
            return sum([linker.n_h_donors for linker in self.linkers])
        except TypeError:
            return None

    def get_num_h_bond_acceptors(self):
        """
        Get the number of hydrogen bond acceptors in a metallocage

        :return: (int)
        """
        try:
            return sum([linker.n_h_acceptors for linker in self.linkers])
        except TypeError:
            return None

    def get_max_escape_sphere(self, basinh=False, max_dist_from_metals=10):
        """
        Get the maximum radius of a sphere that can escape from the centroid of the cage – will iterate through all
        theta/phi

        :param basinh: (bool) Find the true maximum escape sphere by basin hopping on the surface
        :param max_dist_from_metals: (float) Distance in Å on top of the average M-M distance that will be used for the
                                             search for the maximum escape sphere
        :return: (float) Volume of the maximum escape sphere in Å^3
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

        # Get the atom id that the max escape sphere hits into
        sphere_point = spherical_to_cart(r=opt_r, theta=opt_theta_phi[0], phi=opt_theta_phi[1])
        atom_id = np.argmin([np.linalg.norm(coord - sphere_point) for coord in cage_coords])

        radius = max_sphere_escape_r - get_vdw_radii(atom_label=self.xyzs[atom_id][0])
        logger.info(f'Radius of largest sphere that can escape from the cavity = {radius}')

        return (4.0 / 3.0) * np.pi * radius**3

    def _is_linker_reasonable(self, linker):

        if linker is None:
            logger.error(f'Linker was None. Cannot build {self.name}')
            return False

        if linker.xyzs is None or linker.arch is None or linker.name is None:
            logger.error(f'Linker doesn\'t have all the required attributes. Cannot build {self.name}')
            return False

        return True

    def _calc_charge(self):
        logger.info('Calculating the charge on the metallocage')
        self.charge = self.arch.n_metals * self.metal_charge + np.sum(np.array([linker.charge for linker in self.linkers]))
        return None

    def _init_homoleptic_cage(self, linker):
        logger.info(f'Initialising a homoleptic cage')
        self.homoleptic = True

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
        self.heteroleptic = True

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

    def _build(self, with_linker_confs=False, max_per_atom_repulsion=0.001):
        logger.info('Building a cage geometry')
        assert self.homoleptic or self.homoleptic

        xyzs, drs = [], []
        if self.homoleptic:
            linker_conf_list = self.linkers[0].get_ranked_linker_conformers(metal=self.metal)

            # If the linker conformers aren't considered in the build function iterate only over the first (which has
            # the minimum cost function i.e. fits the template best)
            if not with_linker_confs:
                linker_conf_list = linker_conf_list[:1]

            min_repulsion, best_linker = 99999999.9, None
            for linker in linker_conf_list:
                repulsion = 0.0

                for i, template_linker in enumerate(self.cage_template.linkers):

                    linker_xyzs, cost = get_linker_xyzs_to_add_and_cost(linker, template_linker, curr_xyzs=xyzs)
                    repulsion += cost

                    if linker_xyzs is None:
                        logger.error('Failed to get linker')
                        break

                    xyzs += linker_xyzs

                if repulsion < min_repulsion:
                    min_repulsion = repulsion
                    best_linker = deepcopy(linker)

                if repulsion / len(xyzs) < max_per_atom_repulsion:
                    logger.info(f'Total L-L repulsion in building cage is {repulsion:.2f}')
                    self.dr = linker.dr
                    break

                else:
                    xyzs = []

            if len(xyzs) == 0 and best_linker is None:
                logger.error('Could not achieve the required cost threshold for building the cage')
                return

            logger.warning('Failed to reach the threshold. Returning the cage that minimises the L-L repulsion')
            self.dr = best_linker.dr
            xyzs = []
            for i, template_linker in enumerate(self.cage_template.linkers):
                linker_xyzs, _ = get_linker_xyzs_to_add_and_cost(best_linker, template_linker, curr_xyzs=xyzs)
                xyzs += linker_xyzs

        if self.heteroleptic:
            logger.critical('NOT IMPLEMENTED YET')
            exit()

        # Add the metals from the template shifted by dr
        for metal in self.cage_template.metals:
            metal_coord = self.dr * metal.shift_vec / np.linalg.norm(metal.shift_vec) + metal.coord
            xyzs.append([self.metal] + metal_coord.tolist())

        if len(xyzs) != self.arch.n_metals + np.sum(np.array([linker.n_atoms for linker in self.linkers])):
            logger.error('Failed to build a cage')
            return None

        self.set_xyzs(xyzs)
        return None

    def __init__(self, linker=None, metal=None, metal_charge=0, linkers=None, solvent=None, mult=1, name='cage'):
        """
        Metallocage object. Inherits from cgbind.molecule.BaseStruct

        :ivar self.metal: (str)
        :ivar self.linkers: (list(Linker object))
        :ivar self.dr: (float)
        :ivar self.arch: (Arch object)
        :ivar self.cage_template: (Template object)
        :ivar self.m_ids: (list(int))
        :ivar self.metal_charge: (int)

        :param name: (str) Name of the cage
        :param solvent: (str)
        :param linker: (Linker object) Linker to initialise a homoleptic metallocage
        :param linkers: (list(Linker object)) List of Linkers to inialise a metallocage
        :param metal: (str) Atomic symbol of the metal
        :param metal_charge: (int) Formal charge on the metal atom/ion
        :param mult: (int) Total spin multiplicity of the cage
        """
        super(Cage, self).__init__(name=name, charge=0, mult=mult, xyzs=None, solvent=solvent)

        logger.info(f'Initialising a Cage object')

        self.metal = metal
        self.linkers = None
        self.dr = None
        self.arch = None
        self.cage_template = None
        self.m_ids = None
        self.metal_charge = int(metal_charge)

        self.reasonable_geometry = False
        self.homoleptic = False
        self.heteroleptic = False

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
        self._build(with_linker_confs=True)

        if self.xyzs is None:
            self.reasonable_geometry = False
            print_output('Cage build for', self.name, 'Failed')
            return

        self.m_ids = self.get_metal_atom_ids()
        self.reasonable_geometry = is_geom_reasonable(self.xyzs)

        print_output('Cage', self.name, 'Built')
