from cgbind.log import logger
import networkx as nx
import numpy as np
from cgbind.architectures import M2L4
from cgbind.architectures import M4L6
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

        # Order the x_motifs according to the centroid – coord distance: smallest -> largest

        return [sorted(x_motif, key=centroid_atom_distance) for x_motif in x_motifs]

    def _find_centroid(self):

        metal_coords = [xyz2coord(xyz) for xyz in self.xyzs if self.metal in xyz]
        return np.average(metal_coords, axis=0)

    def __init__(self, arch, mol2_filename):

        self.arch = arch
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






if __name__ == '__main__':

    # template = Template(arch=M2L4, mol2_filename='/Users/tom/repos/cgbind/cgbind/lib/EZEVAI.mol2')
    # from cgbind.input_output import xyzs2xyzfile
    # xyzs2xyzfile(xyzs=template.xyzs, basename='template')
    # print(template.x_motifs)

    template2 = Template(arch=M4L6, mol2_filename='/Users/tom/repos/cgbind/cgbind/lib/SALDIV.mol2')
    from cgbind.input_output import xyzs2xyzfile
    xyzs2xyzfile(xyzs=template2.xyzs, basename='template')
    print(template2.x_motifs)
