import os
from copy import deepcopy
import numpy as np
from cgbind.x_motifs import get_shifted_template_x_motif_coords
from cgbind.build import get_template_fitted_coords_and_cost
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
        logger.info('Getting metal_label atom ids with label {}'.format(self.metal))
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
                logger.error('Could not find any metal_label atoms')

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

        self.xyzs, self.energy = opt_geom(self.xyzs, self.name, charge=self.charge, n_cores=n_cores)

    def singlepoint(self, n_cores=1):
        self.energy = singlepoint(self, n_cores)

    def calc_charge(self):
        return self.arch.n_metals * self.metal_charge + self.arch.n_linkers * self.linker.charge

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
        for template_linker in self.cage_template.linkers:

            new_linker = deepcopy(self.linker)
            shifted_coords = get_shifted_template_x_motif_coords(linker_template=template_linker, dr=self.dr)
            x_coords = [new_linker.coords[atom_id] for motif in new_linker.x_motifs for atom_id in motif.atom_ids]

            linker_coords, _ = get_template_fitted_coords_and_cost(linker=self.linker, template_x_coords=shifted_coords,
                                                                   coords_to_fit=x_coords)

            xyzs += [[new_linker.xyzs[i][0]] + linker_coords[i].tolist() for i in range(new_linker.n_atoms)]

        if len(xyzs) != self.arch.n_metals + self.arch.n_linkers * self.linker.n_atoms:
            logger.error('Failed to build a cage')
            return None

        return xyzs

    def __init__(self, linker, metal='Pd', metal_charge=0, name='cage'):
        """
        Initialise a cage object
        :param linker: Linker object
        :param metal: (str)
        :param metal_charge: (int)
        :param name: (str)
        """
        logger.info(f'Initialising a Cage object for {linker.name}')

        self.name = name
        self.metal = metal
        self.linker = linker

        if linker.xyzs is None or linker.arch is None:
            self.reasonable_geometry = False
            logger.error('Linker has no xyzs. Can\'t build a cage')
            return

        self.arch = linker.arch
        self.cage_template = linker.cage_template
        self.dr = linker.dr
        self.metal_charge = metal_charge
        self.charge = self.calc_charge()

        self.reasonable_geometry = True
        self.energy, self.xyzs, self.m_ids = None, None, None
        self.xyzs = self.build()

        if self.xyzs is None:
            self.reasonable_geometry = False
            print_output('Cage build for', self.name, 'Failed')
            return

        self.m_ids = self.get_metal_atom_ids()
        self.reasonable_geometry = is_geom_reasonable(self.xyzs)

        print_output('Cage', self.name, 'Built')
