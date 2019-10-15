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
from cgbind.x_motifs import find_x_motifs
from cgbind.x_motifs import check_x_motifs


def find_mols_in_xyzs(xyzs, allow_same=False):
    """
    From a list of xyzs determine the bonds in the system, thus the distinct molecules in the system
    :param xyzs: (list(list)) standard xyzs
    :param allow_same: (bool) add only the unique molecules (False) or add every molecule (True)
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
        if mol_atom_labels not in unique_mols or allow_same:
            unique_mols.append(mol_atom_labels)
            unique_mol_ids.append(molecule)

    unique_mol_xyzs = [[xyzs[n] for n in mol_ids] for mol_ids in unique_mol_ids]
    logger.info('Found {} molecule(s)'.format(len(unique_mol_xyzs)))

    return unique_mol_xyzs


def get_template(arch_name='m2l4', folder_path=None):
    """
    Get a template defined by the architectures arch_name. The saved objects are in
    :param arch_name: (str)
    :param folder_path: (str) Path to the folder where the templates are stored. Has a default
    :return: (object) Template
    """
    obj_name = arch_name + '.obj'
    logger.info(f'Getting a template object from {obj_name}')

    if folder_path is None:
        folder_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib')

    try:
        with open(os.path.join(folder_path, obj_name), 'rb') as pickled_file:
            return pickle.load(pickled_file)

    except IOError:
        logger.critical('Could not retrieve template')
        exit()


class Linker:

    def __init__(self, xyzs, x_atoms):
        self.xyzs = xyzs
        self.x_atoms = x_atoms
        self.coords = xyz2coord(xyzs)
        self.bonds = get_xyz_bond_list(xyzs=self.xyzs)
        self.centroid = np.average(self.coords, axis=0)

        self.x_motifs = find_x_motifs(self)
        check_x_motifs(self)                                        # check that the x_motifs are the same length -
        self.len_x_motif = len(self.x_motifs[0])                    # only applicable to symmetric structures


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

    def _find_linker_mols_and_x_atoms(self):
        logger.info('Stripping the metals from the structure')
        xyzs_no_metals = [xyz for xyz in self.xyzs if self.metal not in xyz]

        logger.info('Finding the distinct linker molecules ')
        linkers_xyzs = find_mols_in_xyzs(xyzs=xyzs_no_metals, allow_same=True)

        linkers_x_atoms = []
        # Add the x_atoms which are contained within each linker, that were found bonded to each metal
        for xyzs in linkers_xyzs:
            coords = xyz2coord(xyzs)
            linker_x_atoms = []

            for i, coord in enumerate(coords):
                for donor_atom_id in self.x_atoms:
                    if list(coord) == list(self.coords[donor_atom_id]):
                        linker_x_atoms.append(i)
                        break

            linkers_x_atoms.append(linker_x_atoms)

        logger.info(f'Found {len(linkers_xyzs)} linkers each with {len(linker_x_atoms)} donor atoms')
        return linkers_xyzs, linkers_x_atoms

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
        for (i, j) in self.bonds:
            if i in self.metal_atoms:
                donor_atoms.append(j)

            if j in self.metal_atoms:
                donor_atoms.append(i)

        logger.info(f'Found {len(donor_atoms)} donor atoms in the structure')
        return donor_atoms

    def _find_centroid(self):

        metal_coords = [xyz2coord(xyz) for xyz in self.xyzs if self.metal in xyz]
        return np.average(metal_coords, axis=0)

    def save_template(self, folder_path=None):
        logger.info('Saving metallocage template')

        if folder_path is None:
            folder_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib')

        with open(os.path.join(folder_path, self.arch_name + '.obj'), 'wb') as pickled_file:
            pickle.dump(self, file=pickled_file)

        return None

    def __init__(self, arch_name, mol2_filename):

        self.arch_name = arch_name
        self.all_xyzs = mol2file_to_xyzs(filename=mol2_filename)
        self.mols_xyzs = find_mols_in_xyzs(xyzs=self.all_xyzs)
        self.xyzs = None                                            # Set in find_metallocage_mol()
        self.metal = None                                           # Set in find_metallocage_mol(). Limited to 1 metal

        self._find_metallocage_mol()
        self.metal_atoms = self._find_metal_atom_ids()              # metal atom ids
        self.bonds = get_xyz_bond_list(xyzs=self.xyzs)
        self.x_atoms = self._find_donor_atoms()                     # donor atom ids

        self.coords = xyz2coord(xyzs=self.xyzs)
        self.centroid = self._find_centroid()                       # centroid of the cage ~ average metal coordinate

        self.distance_matrix = calc_distance_matrix(xyzs=self.xyzs)

        linkers_xyzs, linkers_x_atoms = self._find_linker_mols_and_x_atoms()
        self.linkers = [Linker(linkers_xyzs[i], linkers_x_atoms[i]) for i in range(len(linkers_xyzs))]
