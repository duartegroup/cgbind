from cgbind.log import logger
import numpy as np
import itertools
from cgbind.config import Config
from multiprocessing import Pool
from cgbind.molecule import Molecule
from cgbind.architectures import archs
from cgbind.atoms import heteroatoms
from cgbind.geom import xyz2coord
from cgbind.build import get_new_linker_and_cost
from cgbind.templates import get_template
from cgbind.x_motifs import find_x_motifs
from cgbind.x_motifs import check_x_motifs
from cgbind.x_motifs import sort_x_motifs


class Linker(Molecule):

    def _find_possible_donor_atoms(self):
        return [i for i in range(self.n_atoms) if self.xyzs[i][0] in heteroatoms]

    def _strip_possible_x_motifs_on_connectivity(self):
        """
        For a list of x motifs remove those which don't have the same number of atoms as the template linkers
        :return: (list) new list of x motifs
        """
        logger.info('Stripping x motifs from liker based on the n atoms in each must = template')
        return [x_motif for x_motif in self.x_motifs
                if x_motif.n_atoms == self.cage_template.linkers[0].x_motifs[0].n_atoms]

    def _set_arch(self, arch_name):

        for arch in archs:
            if arch_name == arch.name:
                self.arch = arch
        return None

    def is_planar(self):
        """
        Determine if the linker is planar

        :return: (bool)
        """

        logger.info('Determining if the linker is planar')

        if len(self.x_motifs) != 2:
            logger.error('Could not calculate planarity of the linker. Only implemented for ')
            return True

        x_coords = [self.coords[atom_id] for motif in self.x_motifs for atom_id in motif.atom_ids]
        x_motifs_centroid = np.average(x_coords, axis=0)

        # Calculate the plane containing the centroid and two atoms in the x motifs
        v1 = self.coords[self.x_motifs[0].atom_ids[0]] - x_motifs_centroid        # First atom in the first x motif
        v2 = self.coords[self.x_motifs[1].atom_ids[-1]] - x_motifs_centroid       # Last atom in the second x motif

        norm = np.cross(v1, v2)
        a, b, c = norm
        d = np.dot(norm, x_motifs_centroid)

        logger.info(f'Equation of the plane is {a}x + {b}y + {c}z + d')

        # Calculate the sum of signed distances, which will be 0 if the linker is flat
        sum_dist = 0
        for coord in self.coords:
            dist = (np.dot(coord, norm) + d) / np.linalg.norm(norm)
            sum_dist += dist

        rel_sum_dist = sum_dist / self.n_atoms
        logger.info(f'Relative sum of distances was {rel_sum_dist}')

        if np.abs(rel_sum_dist) > 1:
            logger.info('Sum of signed plane-coord distances was highly != 0. Linker is probably not planar')
            return False

        else:
            logger.info('Linker is planar')
            return True

    def get_ranked_linker_conformers(self, metal=None):
        """
        For this linker, return a list of Linker objects with appropriate .xyzs, .dr and .x_motifs attributes ordered
        by their cost function low -> high i.e. good to bad. This will loop through all the conformers and the possible
        combinations of x motifs in the linker. Linker.dr controls how large the template needs to be to make
        the best fit

        :param metal: (str) Atomic symbol of the metal
        :return: (list(Linker))
        """
        logger.info('Getting linkers ranked by cost')

        linkers = []

        template_linker = self.cage_template.linkers[0]
        n_x_motifs_in_linker = len(template_linker.x_motifs)

        # For all the possible combinations of x_motifs minimise the SSD between the x_motifs and the template
        # x_motifs. The template needs to be modified to accommodate longer linkers with the same architecture
        x_motifs_list = list(itertools.combinations(self.x_motifs, n_x_motifs_in_linker))

        # Sort the list of x_motifs in the linker by the most favourable M––X interaction
        x_motifs_list = sort_x_motifs(x_motifs_list, linker=self, metal=metal)

        logger.info(f'Have {len(x_motifs_list)*len(self.conf_xyzs)} iterations to do')
        for i, x_motifs in enumerate(x_motifs_list):

            # Execute calculation to get cost of adding a particular conformation to the template in parallel
            logger.info(f'Running with {Config.n_cores} cores. Iteration {i}/{len(x_motifs_list)}')
            with Pool(processes=Config.n_cores) as pool:
                results = [pool.apply_async(get_new_linker_and_cost, (xyzs, self, x_motifs, template_linker))
                           for xyzs in self.conf_xyzs]

                linkers_and_cost_tuples = [res.get(timeout=None) for res in results]

            linkers_and_cost = {}
            for (linker, cost) in linkers_and_cost_tuples:
                linkers_and_cost[linker] = cost

            # Sort this block of linkers the cost function. Not sorted the full list to retain the block structure with
            # X motifs
            linkers += sorted(linkers_and_cost, key=linkers_and_cost.get)

        return linkers

    def __init__(self, arch_name, smiles=None, name='linker', charge=0, n_confs=200, xyzs=None, use_etdg_confs=False):
        """
        Metallocage Linker. Inherits from cgbind.molecule.Molecule

        :param arch_name: (str) Name of the architecture
        :param smiles: (str) SMILES string
        :param name: (str) Linker name
        :param charge: (int)
        :param n_confs: (int) Number of initial conformers to search through
        :param xyzs: (list(list))
        :param use_etdg_confs: (bool) Use a different, sometimes better, conformer generation algorithm
        """

        logger.info('Initialising a Linker object for {}'.format(name))
        initalised_with_xyzs = True if xyzs is not None else False

        super(Linker, self).__init__(smiles=smiles, name=name, charge=charge, n_confs=n_confs, xyzs=xyzs,
                                     use_etdg_confs=use_etdg_confs)

        self.arch = None                                                      #: (Arch object) Metallocage architecture
        self._set_arch(arch_name)

        if self.arch is None:
            logger.error(f'Not a valid architecture. Valid are {[arch.name for arch in archs]}')
            return

        if self.xyzs is None:
            logger.error('Could get xyzs for linker')
            return

        self.cage_template = get_template(arch_name=arch_name)                #: (Template object) Metallocage template

        self.coords = xyz2coord(self.xyzs)                                    #: (list(np.ndarray)) Linker coordinates
        self.centroid = np.average(self.coords, axis=0)                       #: (np.ndarray) Linker centroid ~ COM

        self.x_atoms = self._find_possible_donor_atoms()                      #: (list(int)) List of donor atom ids
        self.x_motifs = find_x_motifs(self, all_possibilities=True)           #: (list(Xmotif object))
        check_x_motifs(self, linker_template=self.cage_template.linkers[0])
        self.x_motifs = self._strip_possible_x_motifs_on_connectivity()
        self.dr = None                                                        #: (float) Template shift distance

        # If the linker has been initialised from xyzs then set conf_xyzs as the xyzs
        if initalised_with_xyzs:
            self.conf_xyzs = [self.xyzs]
