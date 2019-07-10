import os
import numpy as np
from . import m2l4
from . import m4l6
from .log import logger
from .config import Config
from .input_output import xyzfile2xyzs
from .optimisation import opt_geom
from .geom import is_geom_reasonable
from .geom import calc_midpoint
from .geom import xyz2coord


class Cage(object):

    def get_metal_atom_ids(self):
        logger.info('Getting metal atom ids with label {}'.format(self.metal))
        try:
            return [i for i in range(len(self.xyzs)) if self.xyzs[i][0] == self.metal]
        except TypeError or IndexError or AttributeError:
            logger.error('Could not get metal atom ids')
            return

    def get_cavity_vol(self):
        """
        For a cage extract the cavity volume defined as the volume of the largest sphere, centered on the M-M midpoint
        that may be constructed while r < r(midpoint--closest atom)
        :return: Cavity volume in Å^3
        """
        logger.info('Calculating cage cavity vol (will be OVERESTIMATED)')

        cavity_vol = 0.0

        try:

            if self.arch == 'm2l4':
                min_r_midpoint_atom = 999.9
                m_m_midpoint = calc_midpoint(xyz2coord(self.xyzs[self.m_ids[0]]), xyz2coord(self.xyzs[self.m_ids[1]]))
                for i in range(len(self.xyzs)):
                    dist = np.linalg.norm(xyz2coord(self.xyzs[i]) - m_m_midpoint)

                    if dist < min_r_midpoint_atom:
                        min_r_midpoint_atom = dist

                cavity_vol = (4.0 / 3.0) * np.pi * min_r_midpoint_atom**3

            if self.arch == 'm4l6':
                # TODO calc cavity vol for M4L6 cages
                pass

        except TypeError or ValueError or AttributeError:
            logger.error('Could not calculate the cavity volume')
            pass

        return cavity_vol

    def get_m_m_dist(self):
        """
        For a cage calculate the M-M distance
        :return: Distance in Å
        """
        m_m_dist = 0.0

        try:
            if self.arch == 'm2l4':
                logger.info('Calculating the M–M distance in a M2L4 cage')
                m_m_dist = np.linalg.norm(xyz2coord(self.xyzs[self.m_ids[0]]) - xyz2coord(self.xyzs[self.m_ids[1]]))

            if self.arch == 'm4l6':
                logger.critical('Calculating the average M–M distance in a M4L6 cage. NOT IMPLEMENTED')
                exit()

        except TypeError or ValueError or AttributeError:
            logger.error('Could not calculate the M-M distance')
            pass

        return m_m_dist

    def optimise(self, n_cores=1):
        """
        Optimise a cage geometry.
        If there exists an optimised geometry in path_to_opt_struct then cage.xyzs will be set with that geometry.
        :param n_cores: Number of cores to perform the optimisation with (int)
        :return:
        """
        logger.info('Optimising {}'.format(self.name))

        if Config.path_to_opt_struct:
            path_to_opt_geom = os.path.join(Config.path_to_opt_struct, self.name + '.xyz')
            if os.path.exists(path_to_opt_geom):
                logger.info('Found an optimised geometry in path_to_opt_struct')
                self.xyzs = xyzfile2xyzs(xyz_filename=path_to_opt_geom)
        else:
            self.xyzs, self.energy = opt_geom(self.xyzs, self.name, charge=self.charge, n_cores=n_cores)

        return 0

    def optimise_cage_substrate(self, opt_atom_ids=None, n_cores=1):
        """
        Like optimise, but with a cage.substrate species
        :param opt_atom_ids: Atom ids to optimise, possible to just optimise the substrate within the cage,
        providing a considerable computational speedup
        :param n_cores: Number of cores to use for the optimisation
        :return:
        """
        logger.info('Optimising a cage-substrate complex')

        self.cage_substrate_xyzs, self.cage_substrate_energy = opt_geom(self.cage_substrate_xyzs,
                                                                        self.cage_substrate_name,
                                                                        charge=self.charge,
                                                                        opt_atom_ids=opt_atom_ids,
                                                                        n_cores=n_cores)
        return 0

    def add_substrate(self, substrate):
        """
        Add a substrate to a cage.
        The binding mode will be determined by the number of heteroatoms etc. in the case of an M2L4 cage,
        :param substrate: Substrate object
        :return:
        """

        if self.xyzs is None or substrate.xyzs is None:
            logger.error("Cage and/or substrate has no xyzs. Can't add a substrate")
            return 1

        if not Config.suppress_print:
            print("{:<30s}{:<50s}{:>10s}".format('Addition of ', substrate.name, 'Running'))

        self.cage_substrate_name = self.name + '_' + substrate.name
        self.substrate = substrate
        self.cage_substrate_xyzs = self.add_substrate_xyzs(substrate)
        self.charge += substrate.charge
        if self.cage_substrate_xyzs:
            self.cage_subst_reasonable_geometry = is_geom_reasonable(self.cage_substrate_xyzs)
            if not Config.suppress_print:
                print("{:<30s}{:<50s}{:>10s}".format('', substrate.name, 'Done'))
                return 0

        return 1

    def add_substrate_xyzs(self, substrate):
        """
        Add the substrate to the cavity in a mode defined by the number of heteroatoms.

        :param substrate:
        :return:
        """
        logger.info('Adding substrate xyzs to cage')

        xyzs = None

        if self.arch == 'm2l4':
            if substrate.x_x_atom_ids:
                if all(substrate.x_x_atom_ids) is not None:
                    xyzs = m2l4.add_substrate_x_x(self, substrate)
            elif substrate.x_atom_ids:
                xyzs = m2l4.add_substrate_x(self, substrate)
            else:
                xyzs = m2l4.add_substrate_com(self, substrate)

        if self.arch == 'm4l6':
            # TODO add substrate with M4L6 cage
            logger.critical('Adding a substrate to an M4L6 cage. NOT IMPLEMENTED YET')
            exit()

        return xyzs

    def __init__(self, linker, metal='Pd', total_charge=4, name='cage', arch='m2l4'):
        """
        Initialise a cage object
        :param linker: Linker object
        :param metal: Metal atom label (str)
        :param total_charge: Total charge on the cage i.e. metals + linkers (int)
        :param name: Name of the metallocage (str)
        :param arch: Cage architecture. Currently only 'm2l4' or 'm4l4' are supported
        """
        logger.info('Initialising a Cage object for {}'.format(linker.name))

        self.name = name
        self.metal = metal
        self.charge = total_charge
        self.arch = arch.lower()
        self.linker = linker
        self.reasonable_geometry = True
        self.cage_subst_reasonable_geometry = False
        self.energy, self.xyzs, self.m_ids = None, None, None
        self.substrate, self.cage_substrate_name, self.substrate_atom_ids = None, None, None
        self.cage_substrate_xyzs, self.cage_substrate_energy = None, None

        if not linker.xyzs:
            self.reasonable_geometry = False
            logger.error('Linker has no xyzs. Can\'t build a cage')
            return

        if self.arch == 'm2l4':
            self.xyzs = m2l4.build(self, linker)

        elif self.arch == 'm4l6':
            self.xyzs = m4l6.build(self, linker)

        else:
            logger.critical("Couldn't build a cage with architecture {}. NOT IMPLEMENTED YET".format(self.arch))
            exit()

        if self.xyzs is None:
            self.reasonable_geometry = False
            logger.error("Couldn't build a cage")
            return

        self.m_ids = self.get_metal_atom_ids()
        self.reasonable_geometry = is_geom_reasonable(self.xyzs)

        if not Config.suppress_print:
            print("{:<30s}{:<50s}{:>10s}".format('Cage', self.name, 'Built'))
