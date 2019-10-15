from cgbind.log import logger
import networkx as nx
import numpy as np
import pickle
import os
from cgbind.input_output import mol2file_to_xyzs
from cgbind.bonds import get_xyz_bond_list
from cgbind.atoms import metals
from cgbind.geom import xyz2coord
from cgbind.geom import calc_distance_matrix


def find_mols_in_xyzs(xyzs):
    """
    From a list of xyzs determine the bonds in the system, thus the distinct molecules in the system
    :param xyzs: (list(list)) standard xyzs
    :return: (list(xyzs))
    """

    logger.info('Finding the distinct molecules in the system')
    # Get the 'molecules' for which each atom is contained in
    full_graph = nx.Graph()
    [full_graph.add_node(n, atom_label=xyzs[n][0]) for n in range(len(xyzs))]
    bond_list = get_xyz_bond_list(xyzs=xyzs)

    for (u, v) in bond_list:
        full_graph.add_edge(u, v)

    unique_mols, unique_mol_ids = [], []
    connected_molecules = [list(mol) for mol in nx.connected_components(full_graph)]

    for molecule in connected_molecules:
        mol_atom_labels = sorted([xyzs[n][0] for n in molecule])
        if mol_atom_labels not in unique_mols:
            unique_mols.append(mol_atom_labels)
            unique_mol_ids.append(molecule)

    unique_mol_xyzs = [[xyzs[n] for n in mol_ids] for mol_ids in unique_mol_ids]
    logger.info('Found {} distinct molecules'.format(len(unique_mol_xyzs)))

    return unique_mol_xyzs


class Template:

    def _find_metallocage_mol(self):
        """
        From a list of distinct molecules find the metallocage. This is assumed to be the molecule with the highest
        frequency of metal atoms

        :return: None
        """
        mol_metals_and_freqs, metal = [], None

        mols_xyzs = find_mols_in_xyzs(self.all_xyzs)
        for xyzs in mols_xyzs:
            metals_and_freq = dict.fromkeys(metals, 0)

            for atom_label, _, _, _ in xyzs:
                if atom_label in metals:
                    metals_and_freq[atom_label] += 1

            # Add the maximum frequency that any metal arises in the structure
            metal = max(metals_and_freq, key=metals_and_freq.get)
            freq = metals_and_freq[metal]
            mol_metals_and_freqs.append((metal, freq))

            logger.info(f'Max metal frequencies in molecules are {metal} with n = {freq}')

        mol_id_with_max_metals, max_freq, max_metal = 0, 0, None
        for i, (metal, freq) in enumerate(mol_metals_and_freqs):
            if freq > max_freq:
                max_freq = freq
                mol_id_with_max_metals = i
                max_metal = metal

        self.xyzs = self.mols_xyzs[mol_id_with_max_metals]
        self.metal = max_metal
        logger.info(f'Set metal as {self.metal}')
        return None

    def _find_metal_atom_ids(self):
        logger.info('Getting metal atom ids with label {}'.format(self.metal))
        try:
            return [i for i in range(len(self.xyzs)) if self.xyzs[i][0] == self.metal]
        except TypeError or IndexError or AttributeError:
            logger.error('Could not get metal atom ids. Returning None')

        return None

    def _find_donor_atoms(self):
        """
        Find the donor atoms or 'x_atoms' in a metallocage. Will be bonded to the metal with a bond distance up to
        1.1 x the value defined in cgbind.bonds.

        :return: (list(int))
        """
        logger.info('Getting the donor (x) atoms in a structure')

        donor_atoms = []
        for (i, j) in self.bond_list:
            if i in self.metal_atoms:
                donor_atoms.append(j)

            if j in self.metal_atoms:
                donor_atoms.append(i)

        logger.info(f'Found {len(donor_atoms)} donor atoms in the structure')
        return donor_atoms

    def _find_x_motifs(self):
        """
        Find the X motifs in a structure which correspond to the X atoms and their nearest neighbours. These motifs will
        be searched in a linker object and the RMSD minimised

        :return: (list(list(int)))
        """
        def centroid_atom_distance(atom_i):
            return np.linalg.norm(self.coords[atom_i] - self.centroid)

        x_motifs = []

        for donor_atom in self.x_atoms:
            x_motif = []

            # Add all the atoms that are connected to the donor atom
            for (i, j) in self.bond_list:
                if donor_atom == i and j not in self.metal_atoms and j not in x_motif:
                    x_motif.append(j)
                if donor_atom == j and i not in self.metal_atoms and j not in x_motif:
                    x_motif.append(i)

            logger.info(f'X atom {donor_atom} had {len(x_motif)} nearest neighbours')
            x_motif.append(donor_atom)
            x_motifs.append(x_motif)

        # Combine x_motifs that are bonded, thus should be considered a single motif
        bonded_x_motif_sets = []
        for n, x_motif_i in enumerate(x_motifs):
            bonded_x_motif = x_motif_i.copy()

            # Loop over all other motifs that are not the same
            for m, x_motifs_j in enumerate(x_motifs):
                if n != m:

                    for (i, j) in self.bond_list:
                        # Check that there is a bond between x_motif i and x_motif j
                        if (i in bonded_x_motif and j in x_motifs_j) or (j in bonded_x_motif and i in x_motifs_j):

                            bonded_x_motif += x_motifs_j
                            break

            bonded_x_motif_sets.append(set(bonded_x_motif))

        # Add the largest set to the bonded_x_motifs as a list. Some motifs will be missed due to the oder in which
        # they're added

        largest_unique_bonded_x_motif_sets = []
        # Sort the sets according to size, then don't add identical sets or subsets
        for x_motif in reversed(sorted(bonded_x_motif_sets, key=len)):

            unique = True
            for unique_x_motif in largest_unique_bonded_x_motif_sets:

                # If the motif is already in the unique list then don't append, nor if the motif is a subset
                if x_motif == unique_x_motif or x_motif.issubset(unique_x_motif):
                    unique = False
                    break

            if unique:
                largest_unique_bonded_x_motif_sets.append(x_motif)

        logger.info(f'Found {len(largest_unique_bonded_x_motif_sets)} X motifs in the template')

        # Order the x_motifs according to the centroid – coord distance: smallest -> largest
        return [sorted(list(x_motif), key=centroid_atom_distance) for x_motif in largest_unique_bonded_x_motif_sets]

    def _find_centroid(self):

        metal_coords = [xyz2coord(xyz) for xyz in self.xyzs if self.metal in xyz]
        return np.average(metal_coords, axis=0)

    def _check_x_motifs(self):

        if not all([len(motif) == len(self.x_motifs[0]) for motif in self.x_motifs]):
            logger.critical('Found x motifs in the structure that have different number of atoms')
            exit()

        logger.info(f'Number of atoms in the x motifs is {len(self.x_motifs[0])}')
        return None

    def save_template(self, folder_path=None):
        logger.info('Saving metallocage template')

        if folder_path is None:
            folder_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib')

        with open(os.path.join(folder_path, self.name + '.obj'), 'wb') as pickled_file:
            pickle.dump(self, file=pickled_file)

        return None

    def __init__(self, name, mol2_filename):

        self.name = name
        self.all_xyzs = mol2file_to_xyzs(filename=mol2_filename)
        self.mols_xyzs = find_mols_in_xyzs(xyzs=self.all_xyzs)
        self.xyzs = None                                            # Set in find_metallocage_mol()
        self.metal = None                                           # Set in find_metallocage_mol(). Limited to 1 metal

        self._find_metallocage_mol()
        self.metal_atoms = self._find_metal_atom_ids()
        self.bond_list = get_xyz_bond_list(xyzs=self.xyzs)
        self.x_atoms = self._find_donor_atoms()

        self.coords = xyz2coord(xyzs=self.xyzs)
        self.centroid = self._find_centroid()

        self.distance_matrix = calc_distance_matrix(xyzs=self.xyzs)
        self.x_motifs = self._find_x_motifs()
        self._check_x_motifs()
        self.len_x_motif = len(self.x_motifs[0])


if __name__ == '__main__':

    template = Template(name='m2l4', mol2_filename='lib/EZEVAI.mol2')
    template.save_template()

    template2 = Template(name='m4l6', mol2_filename='lib/GARWUR.mol2')
    template2.save_template()

    template3 = Template(name='m4l6_2', mol2_filename='lib/SALDIV.mol2')
    template3.save_template()
