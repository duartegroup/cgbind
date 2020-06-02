import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import basinhopping, minimize
from cgbind.molecule import BaseStruct
from cgbind.calculations import get_charges
from cgbind.build import build_homoleptic_cage
from cgbind.build import build_heteroleptic_cage
from cgbind.log import logger
from cgbind.atoms import get_vdw_radii
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

        metal_coords = np.array([atom.coord for atom in self.atoms if self.metal == atom.label])
        return np.average(metal_coords, axis=0)

    def get_esp_cube(self, return_min_max=False):
        """
        Get the electrostatic potential (ESP) in a Gaussian .cube format by calculating partial atomic charges using
        XTB (tested with v. 6.2). Calls self.get_charges() and depends on self.xyzs

        :param return_min_max: (bool) Return the minimum and maximum of the ESP along with the cube file lines
                                      evaluated roughly on the VdW surface
        :return: (list) .cube file lines
        """
        esp_lines, (min_esp, max_esp) = get_esp_cube_lines(charges=self.get_charges(), atoms=self.atoms)

        if return_min_max:
            return esp_lines, min_esp, max_esp

        else:
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

    def get_charges(self, estimate=False):
        """
        Get the partial atomic charges on the cage either using XTB or estimate using no polarisation i.e. the metals
        retain their full charge and the linker charges are estimated using the Gasteiger scheme in RDKit

        :param estimate: (bool)
        :param guess: (bool) Guess the charges based on the electronegativity
        :return: (function) calculations.get_charges(self)
        """
        if estimate:

            charges = []
            for linker in self.linkers:
                linker_charges = linker.get_charges(estimate=estimate)
                charges += linker_charges

            # Metals are added last
            charges += [self.metal_charge] * self.arch.n_metals

            return charges

        return get_charges(self)

    def get_metal_atom_ids(self):
        """
        Get the atom ids of the self.metal atoms in the xyzs

        :return: (list(int))
        """

        logger.info(f'Getting metal_label atom ids with label {self.metal}')
        if self.n_atoms == 0:
            logger.error('Could not get metal atom ids. xyzs were None')
            return None

        try:
            return [i for i in range(self.n_atoms) if self.atoms[i].label == self.metal]
        except TypeError or IndexError or AttributeError:
            logger.error('Could not get metal label atom ids. Returning None')
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
            for i in range(self.n_atoms):
                dist = np.linalg.norm(self.atoms[i].coord - centroid)
                if dist < min_centriod_atom_dist:
                    min_centriod_atom_dist = dist
                    min_atom_dist_id = i

        except TypeError or ValueError or AttributeError:
            pass

        if min_atom_dist_id is not None:
            vdv_radii = get_vdw_radii(atom=self.atoms[min_atom_dist_id])
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
                        dist = np.linalg.norm(self.atoms[self.m_ids[m_id_i]].coord -
                                              self.atoms[self.m_ids[m_id_j]].coord)
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
            print(self.linkers[0].x_atoms)
            n_donor_atoms = sum([len(linker.x_atoms) for linker in self.linkers])

            return max(sum([linker.n_h_acceptors for linker in self.linkers]) - n_donor_atoms, 0)
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

        centroid = self.get_centroid()
        cage_coords = self.get_coords()
        cage_coords = np.array([coord - centroid for coord in cage_coords])

        # For a distance from the origin (the cage centroid) calculate the largest sphere possible without hitting atoms
        opt_theta_phi, opt_r = np.zeros(2), 0.0
        for r in np.linspace(0.0, avg_m_m_dist + max_dist_from_metals, 30):
            if basinh:
                opt = basinhopping(get_max_sphere_negative_radius, x0=opt_theta_phi, stepsize=1.0, niter=5,
                                   minimizer_kwargs={'args': (r, cage_coords), 'method': 'BFGS'})
            else:
                opt = minimize(get_max_sphere_negative_radius, x0=opt_theta_phi, args=(r, cage_coords), method='BFGS')

            opt_theta_phi = opt.x

            # This is the correct way round because we want the largest sphere that CAN escape
            if -opt.fun < max_sphere_escape_r:
                max_sphere_escape_r = -opt.fun
                opt_r = r

        # Get the atom id that the max escape sphere hits into
        sphere_point = spherical_to_cart(r=opt_r, theta=opt_theta_phi[0], phi=opt_theta_phi[1])
        atom_id = np.argmin([np.linalg.norm(coord - sphere_point) for coord in cage_coords])

        radius = max_sphere_escape_r - get_vdw_radii(atom=self.atoms[atom_id])
        logger.info(f'Radius of largest sphere that can escape from the cavity = {radius}')

        return (4.0 / 3.0) * np.pi * radius**3

    def _is_linker_reasonable(self, linker):

        if linker is None:
            logger.error(f'Linker was None. Cannot build {self.name}')
            return False

        if linker.n_atoms == 0 or linker.arch is None or linker.name is None:
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

        if self.name == 'cage':
            # Only override the default name
            self.name = 'cage_' + linker.name

        self.arch = linker.arch
        self.linkers = [linker for _ in range(linker.arch.n_linkers)]
        self.cage_template = linker.cage_template

        return None

    def _init_heteroleptic_cage(self, linkers):
        logger.info(f'Initialising a heteroleptic cage')
        self.heteroleptic = True

        if not all([self._is_linker_reasonable(linker) for linker in linkers]):
            logger.error('Not all linkers were reasonable')
            return

        if not all([linker.arch.name == linkers[0].arch.name for linker in linkers]):
            logger.error('Linkers had different architectures, not building a cage')
            return

        if self.name == 'cage':
            # Only override the default name
            self.name = 'cage_' + '_'.join([linker.name for linker in linkers])

        self.arch = linkers[0].arch
        self.linkers = linkers
        self.cage_template = linkers[0].cage_template

        return None

    def _build(self, max_cost):
        logger.info('Building a cage geometry')
        assert self.homoleptic or self.heteroleptic

        if self.homoleptic:
            build_homoleptic_cage(self, max_cost)

        if self.heteroleptic:
            build_heteroleptic_cage(self, max_cost)

        if self.reasonable_geometry:
            if self.n_atoms != self.arch.n_metals + np.sum(np.array([linker.n_atoms for linker in self.linkers])):
                logger.error('Failed to build a cage')
                self.reasonable_geometry = False
                return None

        return None

    def __init__(self, linker=None, metal='M', metal_charge=0, linkers=None, solvent=None, mult=1, name='cage', max_cost=5):
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
        :param max_cost: (float) Acceptable ligand-ligand repulsion to accommodate in metallocage construction
        """
        super(Cage, self).__init__(name=name, charge=0, mult=mult, filename=None, solvent=solvent)

        logger.info(f'Initialising a Cage object')

        self.metal = str(metal)
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

        self.reasonable_geometry = False
        self._build(max_cost=max_cost)

        self.m_ids = self.get_metal_atom_ids()
        logger.info(f'Generated cage successfully. Geometry is reasonable: {self.reasonable_geometry}')
