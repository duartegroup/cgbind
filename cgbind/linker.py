from cgbind.log import logger
import numpy as np
import itertools
from cgbind.config import Config
from cgbind.molecule import Molecule
from cgbind.architectures import archs
from cgbind.exceptions import *
from cgbind.atoms import heteroatoms
from cgbind.atoms import get_max_valency
from cgbind.build import calc_conformer_cost
from cgbind.templates import get_template
from cgbind.x_motifs import find_x_motifs
from cgbind.x_motifs import check_x_motifs
from cgbind.x_motifs import sort_x_motifs


class Linker(Molecule):

    def __eq__(self, other):
        """
        Linkers are the same if their SMILES strings are identical, otherwise very close centroids should serve as a
        unique measure for linkers with identical (up to numerical precision) coordinates

        :param other: (Linker)
        :return:
        """

        if self.smiles is not None:
            return self.smiles == other.smiles

        dist = np.linalg.norm(self.com - other.centroid)
        return True if dist < 1E-9 else False

    def __hash__(self):
        return hash((self.name, self.smiles, self.charge, self.com[0]))

    def _find_possible_donor_atoms(self):
        """
        For the atoms in the linker find all those capable of donating a 'lone pair' to a metal i.e. being a donor/
        'X-atom'

        :return: (list(int)) donor atoms
        """
        donor_atom_ids = []

        for atom_id, atom in enumerate(self.atoms):

            if atom.label in heteroatoms:
                max_valency = get_max_valency(atom=atom)
                n_bonds = 0

                for bond in self.bonds:
                    if atom_id in bond:
                        n_bonds += 1

                # If the number of bonds is lower than the max valancy for that atom then there should be a lone pair
                # and a donor atom has been found
                if n_bonds < max_valency:
                    donor_atom_ids.append(atom_id)

        logger.info(f'Found {len(donor_atom_ids)} possible X atoms')
        return donor_atom_ids

    def _strip_possible_x_motifs_on_connectivity(self):
        """
        For a list of x motifs remove those which don't have the same number of atoms as the template linkers
        :return: (list) new list of x motifs
        """
        logger.info('Stripping x motifs from liker based on the n atoms in each must = template')
        self.x_motifs = [x_motif for x_motif in self.x_motifs
                         if x_motif.n_atoms == self.cage_template.linkers[0].x_motifs[0].n_atoms]

        logger.info(f'Current number of donor atoms = {len(self.x_motifs)}')
        self.x_atoms = [atom_id for x_motif in self.x_motifs for atom_id in x_motif.atom_ids if atom_id in self.x_atoms]
        logger.info(f'New number of donor atoms = {len(self.x_motifs)}')
        return None

    def _set_arch(self, arch_name):

        for arch in archs:
            if arch_name.lower() == arch.name.lower():
                self.arch = arch

        if self.arch is None:
            raise ArchitectureNotFound(message=f'\nAvailable architecture names are {[arch.name for arch in archs]}')

        return None

    def _check_structure(self):
        if self.n_atoms == 0:
            logger.error('Could get atoms for linker')
            raise NoXYZs

        return None

    def set_ranked_linker_conformers(self, metal=None, n=0):
        """
        For this linker, return a list of Linker objects with appropriate .xyzs, .dr and .x_motifs attributes ordered
        by their cost function low -> high i.e. good to bad. This will loop through all the conformers and the possible
        combinations of x motifs in the linker. Linker.dr controls how large the template needs to be to make
        the best fit

        :param metal: (str) Atomic symbol of the metal
        :param n: (int) linker number in the template
        :return: (list(Linker))
        """
        logger.info('Getting linkers ranked by cost')

        template_linker = self.cage_template.linkers[n]
        n_x_motifs_in_linker = len(template_linker.x_motifs)

        # For all the possible combinations of x_motifs minimise the SSD between the x_motifs and the template
        # x_motifs. The template needs to be modified to accommodate longer linkers with the same architecture
        x_motifs_list = list(itertools.combinations(self.x_motifs, n_x_motifs_in_linker))

        # Sort the list of x_motifs in the linker by the most favourable M––X interaction
        x_motifs_list = sort_x_motifs(x_motifs_list, linker=self, metal=metal)

        logger.info(f'Have {len(x_motifs_list)*len(self.conformers)} iterations to do')

        conformers = []
        for i, x_motifs in enumerate(x_motifs_list):

            # Execute calculation to get cost of adding a particular conformation to the template in parallel
            logger.info(f'Running with {Config.n_cores} cores. Iteration {i+1}/{len(x_motifs_list)}')
            logger.disabled = True

            # Sort this block of linkers the cost function. Not sorted the full list to retain the block structure with
            # X motifs
            chunk = [calc_conformer_cost(conf, self, x_motifs, template_linker) for conf in self.conformers]
            conformers += sorted(chunk, key=lambda conf: conf.cost)

            # Renable the logging that would otherwise dominate with lots of conformers and/or Xmotifs
            logger.disabled = False

        # Reset the conformers as those with both different Xmotifs and geometry i.e.
        # len(conformers) = len(x_motifs_list)*len(self.conformers)
        self.conformers = conformers

        return None

    def __init__(self, arch_name, smiles=None, name='linker', charge=0, n_confs=300, filename=None, use_etdg_confs=False):
        """
        Metallocage Linker. Inherits from cgbind.molecule.Molecule

        :param arch_name: (str) Name of the architecture
        :param smiles: (str) SMILES string
        :param name: (str) Linker name
        :param charge: (int)
        :param n_confs: (int) Number of initial conformers to search through
        :param filename: (str)
        :param use_etdg_confs: (bool) Use a different, sometimes better, conformer generation algorithm
        """
        logger.info(f'Initialising a Linker object for {name}')

        self.arch = None                                                      #: (Arch object) Metallocage architecture
        self._set_arch(arch_name)

        # May exit here if the specified architecture is not found
        super(Linker, self).__init__(smiles=smiles, name=name, charge=charge, n_confs=n_confs, filename=filename,
                                     use_etdg_confs=use_etdg_confs)

        self._check_structure()
        self.cage_template = get_template(arch_name=arch_name)                #: (Template object) Metallocage template

        self.x_atoms = self._find_possible_donor_atoms()                      #: (list(int)) List of donor atom ids
        self.x_motifs = find_x_motifs(self)                                   #: (list(Xmotif object))
        check_x_motifs(self, linker_template=self.cage_template.linkers[0])
        self._strip_possible_x_motifs_on_connectivity()
        self.dr = None                                                        #: (float) Template shift distance
