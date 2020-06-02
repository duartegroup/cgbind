from cgbind.log import logger
import networkx as nx
import numpy as np
import pickle
import os
from cgbind.input_output import mol2file_to_atoms
from cgbind.bonds import get_bond_list_from_atoms
from cgbind.geom import calc_com
from cgbind.atoms import metals
from cgbind.x_motifs import find_x_motifs
from cgbind.x_motifs import get_maximally_connected_x_motifs
from cgbind.x_motifs import check_x_motifs


def find_mols_in_xyzs(atoms, allow_same=False):
    """
    From a list of xyzs determine the bonds in the system, thus the distinct molecules in the system

    :param atoms: (list(cgbind.atoms.Atom))
    :param allow_same: (bool) add only the unique molecules (False) or add every molecule (True)
    :return: (list(xyzs))
    """

    logger.info('Finding the distinct molecules in the system')
    # Get the 'molecules' for which each atom is contained in
    full_graph = nx.Graph()
    [full_graph.add_node(n, atom_label=atoms[n].label) for n in range(len(atoms))]
    bond_list = get_bond_list_from_atoms(atoms)

    for (u, v) in bond_list:
        full_graph.add_edge(u, v)

    unique_mols, unique_mol_ids = [], []
    connected_molecules = [list(mol) for mol in nx.connected_components(full_graph)]

    for molecule in connected_molecules:
        mol_atom_labels = sorted([atoms[n].label for n in molecule])
        if mol_atom_labels not in unique_mols or allow_same:
            unique_mols.append(mol_atom_labels)
            unique_mol_ids.append(molecule)

    unique_mol_atoms = [[atoms[n] for n in mol_ids] for mol_ids in unique_mol_ids]
    logger.info(f'Found {len(unique_mol_atoms)} molecule(s)')

    return unique_mol_atoms


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


class Metal:

    def __init__(self, label, atom_id, coord):
        """
        Template metal

        :param label: (str) Atomic symbol of the metal
        :param atom_id: (int) Atom id in the Template
        :param coord: (np.ndarray) Coordinate of the atom
        """

        self.label = label                                #: (str) Atomic symbol
        self.atom_id = atom_id                            #: (int)
        self.coord = coord                                #: (np.ndarray) Coordinate of the atom/ion
        self.shift_vec = None                             #: (np.ndarray) Vector to shift by when expanding/contracting


class Linker:

    def __init__(self, atoms, x_atoms):
        """
        Make a template linker object from the corresponding xyzs and the donor atoms which were bonded to the metals
        which form the basis of the cage

        :param atoms: (list(list))
        :param x_atoms: (list(int)) Donor atom ids in the xyzs
        """
        logger.info('Generating a template linker...')
        self.atoms = atoms                                      #: (list(list))
        self.x_atoms = x_atoms                                  #: (list(int)) List of donor atoms in the linker
        self.coords = np.array([atom.coord for atom in atoms])  #: (list(np.ndarray)) Linker coordinates
        self.bonds = get_bond_list_from_atoms(self.atoms)       #: (list(Atom))
        self.com = calc_com(atoms)                              #: (np.ndarray)

        self.x_motifs = find_x_motifs(self)                    #: (list(Xmotif objects)
        self.x_motifs = get_maximally_connected_x_motifs(self.x_motifs, x_atoms=x_atoms)

        check_x_motifs(linker_template=self)            # check that the x_motifs are the same length


class Template:

    def _find_metallocage_mol(self):
        """
        From a list of distinct molecules find the metallocage. This is assumed to be the molecule with the highest
        frequency of metal_label atoms

        :return: atoms, metal_label label
        """
        mol_metals_and_freqs, metal = [], None

        for atoms in self.mols_atoms:
            metals_and_freq = dict.fromkeys(metals, 0)

            for atom in atoms:
                if atom.label in metals:
                    metals_and_freq[atom.label] += 1

            # Add the maximum frequency that any metal_label arises in the structure
            metal = max(metals_and_freq, key=metals_and_freq.get)
            freq = metals_and_freq[metal]
            mol_metals_and_freqs.append((metal, freq))

            logger.info(f'Max metal_label frequencies in molecules are {metal} with n = {freq}')

        mol_id_with_max_metals, max_freq, max_metal = 0, 0, None
        for i, (metal, freq) in enumerate(mol_metals_and_freqs):
            if freq > max_freq:
                max_freq = freq
                mol_id_with_max_metals = i
                max_metal = metal

        logger.info(f'Found metal_label {max_metal}')
        return self.mols_atoms[mol_id_with_max_metals], max_metal

    def _find_linkers(self):
        logger.info('Stripping the metals from the structure')
        atoms_no_metals = [atom for atom in self.atoms if self.metal_label != atom.label]

        logger.info('Finding the distinct linker molecules ')
        linkers_atoms = find_mols_in_xyzs(atoms=atoms_no_metals, allow_same=True)

        linkers = []
        # Add the x_atoms which are contained within each linker, that were found bonded to each metal_label
        for atoms in linkers_atoms:
            coords = np.array([atom.coord for atom in atoms])
            linker_x_atoms = []

            # Iterate through the coordinates until one matches that of the full template
            for i, coord in enumerate(coords):
                for donor_atom_id in self.x_atoms:
                    if list(coord) == list(self.coords[donor_atom_id]):
                        linker_x_atoms.append(i)
                        break

            linkers.append(Linker(atoms=atoms, x_atoms=linker_x_atoms))
            logger.info(f'Linker has {len(linker_x_atoms)} donor atoms')

        logger.info(f'Found {len(linkers_atoms)} linkers each with {len(linker_x_atoms)} donor atoms')
        return linkers

    def _find_metals(self):
        logger.info(f'Getting metals with label {self.metal_label}')
        metals = []
        try:
            for i in range(len(self.atoms)):
                if self.atoms[i].label == self.metal_label:
                    metals.append(Metal(label=self.metal_label, atom_id=i, coord=self.atoms[i].coord))

            return metals

        except (TypeError, IndexError, AttributeError):
            logger.error('Could not get metal_label atom ids. Returning None')

        return None

    def _find_donor_atoms(self):
        """
        Find the donor atoms or 'x_atoms' in a metallocage. Will be bonded to the metal_label with a bond distance up to
        1.1 x the value defined in cgbind.bonds.

        :return: (list(int))
        """
        logger.info('Getting the donor (x) atoms in a structure')

        donor_atoms = []
        for (i, j) in self.bonds:
            if i in [metal.atom_id for metal in self.metals]:
                donor_atoms.append(j)

            if j in [metal.atom_id for metal in self.metals]:
                donor_atoms.append(i)

        logger.info(f'Found {len(donor_atoms)} donor atoms in the structure')
        return donor_atoms

    def _find_centroid(self):

        metal_coords = [atom.coord for atom in self.atoms if self.metal_label == atom.label]
        return np.average(metal_coords, axis=0)

    def _set_shift_vectors(self):
        """
        For the linkers set the shift vectors for the x motifs. These are the centroid -> closest metal_label atom
        vector

        :return: None
        """
        metal_atom_ids = [metal.atom_id for metal in self.metals]

        for linker in self.linkers:
            for x_motif in linker.x_motifs:
                motif_atom_id = x_motif.atom_ids[0]

                closest_metal_id = None
                closest_dist = 99999.9
                for metal_id in metal_atom_ids:
                    dist = np.linalg.norm(linker.coords[motif_atom_id] - self.coords[metal_id])
                    if dist < closest_dist:
                        closest_dist = dist
                        closest_metal_id = metal_id

                shift_vec = self.coords[closest_metal_id] - self.centroid

                # Set the shift vector for the metal closest to this x_motif
                for metal in self.metals:
                    if metal.atom_id == closest_metal_id:
                        metal.shift_vec = shift_vec

                # Set the shift vector for the x motif
                x_motif.shift_vec = shift_vec
                x_motif.r = np.linalg.norm(x_motif.shift_vec)
                x_motif.norm_shift_vec = x_motif.shift_vec / x_motif.r

        return None

    def save_template(self):
        """
        Save the template to ./lib/self.arch_name.obj

        :return: None
        """
        logger.info('Saving metallocage template')

        # Templates will be saved to here/lib/
        folder_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib')

        with open(os.path.join(folder_path, self.arch_name + '.obj'), 'wb') as pickled_file:
            pickle.dump(self, file=pickled_file)

        return None

    def __init__(self, arch_name, mol2_filename):
        """
        Initialise a Template object

        :param arch_name: (str) Name of the architecture
        :param mol2_filename: (str) Name of the .mol2 file to read a structure from e.g. downloaded from the CCDC
        """

        self.arch_name = arch_name

        all_atoms = mol2file_to_atoms(filename=mol2_filename)
        self.mols_atoms = find_mols_in_xyzs(atoms=all_atoms)
        self.atoms, self.metal_label = self._find_metallocage_mol()

        self.metals = self._find_metals()
        self.n_metals = len(self.metals)
        assert self.n_metals > 0
        logger.info(f'Found {self.n_metals} metals')

        self.bonds = get_bond_list_from_atoms(self.atoms)
        self.x_atoms = self._find_donor_atoms()                # donor atom ids

        self.coords = np.array([atom.coord for atom in self.atoms])
        self.centroid = self._find_centroid()                  # centroid of the cage ~ average metal_label coordinate

        self.linkers = self._find_linkers()
        self.n_linkers = len(self.linkers)
        logger.info(f'Found {self.n_linkers} linkers')

        self._set_shift_vectors()
