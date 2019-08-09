import os
import numpy as np
from . import m2l4
from . import m4l6
from .log import logger
from .config import Config
from .input_output import xyzfile2xyzs
from .input_output import print_output
from .optimisation import opt_geom
from .single_point import singlepoint
from .geom import is_geom_reasonable
from .geom import calc_midpoint
from .geom import xyz2coord
from .architectures import M2L4
from .architectures import M4L6


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

            if self.arch == M2L4:
                min_r_midpoint_atom = 999.9
                m_m_midpoint = calc_midpoint(xyz2coord(self.xyzs[self.m_ids[0]]), xyz2coord(self.xyzs[self.m_ids[1]]))
                for i in range(len(self.xyzs)):
                    dist = np.linalg.norm(xyz2coord(self.xyzs[i]) - m_m_midpoint)

                    if dist < min_r_midpoint_atom:
                        min_r_midpoint_atom = dist

                cavity_vol = (4.0 / 3.0) * np.pi * min_r_midpoint_atom**3

            if self.arch == M4L6:
                logger.critical('Calculating the cavity volume in a M4L6 cage. NOT IMPLEMENTED')
                exit()

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
            if self.arch == M2L4:
                logger.info('Calculating the M–M distance in a M2L4 cage')
                m_m_dist = np.linalg.norm(xyz2coord(self.xyzs[self.m_ids[0]]) - xyz2coord(self.xyzs[self.m_ids[1]]))

            if self.arch == M4L6:
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

    def singlepoint(self, n_cores=1):
        self.energy = singlepoint(self, n_cores)

    def calc_charge(self, metal_charge):
        return self.arch.n_metals * metal_charge + self.arch.n_linkers * self.linker.charge

    def __init__(self, linker, metal='Pd', metal_charge=0, name='cage', arch=M2L4):
        """
        Initialise a cage object
        :param linker: Linker object
        :param metal: Metal atom label (str)
        :param metal_charge: Total charge on the cage i.e. metals + linkers (int)
        :param name: Name of the metallocage (str)
        :param arch: Cage architecture. Currently only 'm2l4' or 'm4l4' are supported
        """
        logger.info('Initialising a Cage object for {}'.format(linker.name))

        self.name = name
        self.metal = metal
        self.arch = arch
        self.linker = linker
        self.metal_charge = metal_charge
        self.charge = self.calc_charge(metal_charge)

        self.reasonable_geometry = True
        self.energy, self.xyzs, self.m_ids = None, None, None

        if not linker.xyzs:
            self.reasonable_geometry = False
            logger.error('Linker has no xyzs. Can\'t build a cage')
            return

        if self.arch == M2L4:
            self.xyzs = m2l4.build(self, linker)

        elif self.arch == M4L6:
            self.xyzs = m4l6.build(self, linker)

        else:
            logger.critical("Couldn't build a cage with architecture {}. NOT IMPLEMENTED YET".format(self.arch))
            exit()

        if self.xyzs is None:
            self.reasonable_geometry = False
            print_output('Cage build for', self.name, 'Failed')
            return

        self.m_ids = self.get_metal_atom_ids()
        self.reasonable_geometry = is_geom_reasonable(self.xyzs)

        print_output('Cage', self.name, 'Built')
