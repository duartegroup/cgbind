import os
import numpy as np
from cgbind import m2l4
from cgbind import m4l6
from cgbind.log import logger
from cgbind.config import Config
from cgbind.input_output import xyzfile2xyzs
from cgbind.input_output import print_output
from cgbind.optimisation import opt_geom
from cgbind.single_point import singlepoint
from cgbind.atoms import get_vdw_radii
from cgbind.geom import is_geom_reasonable
from cgbind.geom import xyz2coord



class Cage(object):

    def get_metal_atom_ids(self):
        logger.info('Getting metal atom ids with label {}'.format(self.metal))
        try:
            return [i for i in range(len(self.xyzs)) if self.xyzs[i][0] == self.metal]
        except TypeError or IndexError or AttributeError:
            logger.error('Could not get metal atom ids. Returning None')
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
            metal_coords = [xyz2coord(xyz) for xyz in self.xyzs if self.metal in xyz]
            centroid = np.average(metal_coords, axis=0)

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
                logger.error('Could not find any metal atoms')

        except TypeError or ValueError or AttributeError:
            logger.error('Could not calculate the M-M distance. Returning 0.0')

        return 0.0

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
                self.xyzs = xyzfile2xyzs(filename=path_to_opt_geom)
        else:
            self.xyzs, self.energy = opt_geom(self.xyzs, self.name, charge=self.charge, n_cores=n_cores)

    def singlepoint(self, n_cores=1):
        self.energy = singlepoint(self, n_cores)

    def calc_charge(self, metal_charge):
        return self.arch.n_metals * metal_charge + self.arch.n_linkers * self.linker.charge

    def __init__(self, linker, metal='Pd', metal_charge=0, name='cage'):
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
        self.linker = linker
        self.metal_charge = metal_charge
        self.charge = self.calc_charge(metal_charge)

        self.reasonable_geometry = True
        self.energy, self.xyzs, self.m_ids = None, None, None

        if not linker.xyzs:
            self.reasonable_geometry = False
            logger.error('Linker has no xyzs. Can\'t build a cage')
            return

        # m2l4.build(self, linker)





        if self.xyzs is None:
            self.reasonable_geometry = False
            print_output('Cage build for', self.name, 'Failed')
            return

        self.m_ids = self.get_metal_atom_ids()
        self.reasonable_geometry = is_geom_reasonable(self.xyzs)

        print_output('Cage', self.name, 'Built')
