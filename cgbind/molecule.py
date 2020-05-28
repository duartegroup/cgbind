import numpy as np
from copy import deepcopy
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdPartialCharges
from cgbind.log import logger
from cgbind.exceptions import CgbindCritical
from cgbind.atoms import Atom
from cgbind.geom import is_geom_reasonable
from cgbind.config import Config
from cgbind.bonds import get_bond_list_from_atoms
from cgbind.input_output import get_atoms_from_file
from cgbind.input_output import atoms_to_xyz_file
from cgbind.geom import calc_com
from cgbind import calculations


def extract_conformers_from_rdkit_mol_object(mol_obj, conf_ids):
    """
    Generate xyz lists for all the conformers in conf_ids
    :param mol_obj: Molecule object
    :param conf_ids: (list) list of conformer ids to convert to xyz
    :return: (list(list(cgbind.atoms.Atom)))
    """
    conformers = []

    for i in range(len(conf_ids)):
        mol_block_lines = Chem.MolToMolBlock(mol_obj, confId=conf_ids[i]).split('\n')
        atoms = []

        for line in mol_block_lines:
            split_line = line.split()
            if len(split_line) == 16:
                atom_label, x, y, z = split_line[3], split_line[0], split_line[1], split_line[2]
                atoms.append(Atom(atom_label, float(x), float(y), float(z)))

        conformer = BaseStruct()
        conformer.set_atoms(atoms)
        conformers.append(conformer)

    if len(conformers) == 0:
        logger.critical('Length of conformer xyz list was 0')
        exit()

    return conformers


class BaseStruct:

    def print_xyz_file(self, filename=None):
        """
        Print a .xyz file from self.xyzs provided self.reasonable_geometry is True

        :param filename: (str) Override the default filename
        :return: None
        """
        if not self.reasonable_geometry:
            logger.error('Geometry is not reasonable')

        filename = filename if filename is not None else f'{self.name}.xyz'

        atoms_to_xyz_file(atoms=self.atoms, filename=filename)
        return None

    def get_coords(self):
        return np.array([atom.coord for atom in self.atoms])

    def centre(self):
        """Centre the species at the centre of mass (i.e. so COM = origin)"""

        for atom in self.atoms:
            atom.coord -= self.com

        return None

    def singlepoint(self, method, keywords=None, n_cores=None):
        """
        Perform a single-point energy evaluation using an electronic structure theory method e.g. XTB, ORCA, G09

        :param method: (autode.ElectronicStructureMethod)
        :param keywords: (list(str)) Keywords to use for the ESM e.g. ['SP', 'PBE', 'def2-SVP']
        :param n_cores: (int) Number of cores for the calculation to use
        :return: None
        """
        n_cores = n_cores if n_cores is not None else Config.n_cores
        return calculations.singlepoint(self, method, keywords, n_cores)

    def optimise(self, method, keywords=None, n_cores=1, cartesian_constraints=None):
        """
        Perform a single-point energy evaluation using an electronic structure theory method e.g. XTB, ORCA, G09

        :param method: (autode.ElectronicStructureMethod)
        :param keywords: (list(str)) Keywords to use for the ESM e.g. ['Opt', 'PBE', 'def2-SVP']
        :param n_cores: (int) Number of cores for the calculation to use
        :param cartesian_constraints: (list(int)) List of atom ids to constrain to their current coordinates
        :return: None
        """
        n_cores = n_cores if n_cores is not None else Config.n_cores
        return calculations.optimise(self, method, keywords, n_cores, cartesian_constraints)

    def set_atoms(self, atoms=None, coords=None):
        """
        Set the xyzs of a molecular structure

        :param atoms: (list(cgbind.atoms.Atom))
        :param coords: (np.ndarray) n_atoms x 3 positions of the atoms
        :return: None
        """
        # Reset the atoms in this species using an array of coordinates
        if coords is not None:
            assert type(coords) == np.ndarray
            assert coords.shape == (self.n_atoms, 3)

            # Set the coordinates on a copy of the atoms
            atoms = deepcopy(self.atoms)
            for i, coord in enumerate(coords):
                atoms[i].coord = coord

        # Reset the atoms, number of atoms and the centre of mass
        if atoms is not None:
            assert type(atoms) == list
            assert len(atoms) > 0
            assert hasattr(atoms[0], 'label')
            assert hasattr(atoms[0], 'coord')
            assert len(atoms[0].coord) == 3

            self.atoms = atoms
            self.n_atoms = len(atoms)
            self.com = calc_com(atoms=self.atoms)

            self.reasonable_geometry = is_geom_reasonable(self)
            logger.info(f'Geometry is reasonable: {self.reasonable_geometry}')

        else:
            self.reasonable_geometry = False
            logger.warning('xyzs were None -> n_atoms also None & geometry is *not* reasonable')

        return None

    def __init__(self, name='molecule', charge=0, mult=1, filename=None, solvent=None):
        """
        Base structure class

        :param name: (str)
        :param charge: (int)
        :param mult: (int) Spin multiplicity
        :param filename: (str)
        :param solvent: (str)
        """

        self.name = str(name)                                               #: (str) Name of the structure

        self.n_atoms = 0                                                    #: (int) Number of atoms
        self.atoms = []                                                     #: (list(list)) Geometry of the structure
        self.com = None                                                     #: (np.ndarray) Center of mass (x, y, z)

        # If initialised from a structure file then set the atoms
        if filename is not None:
            self.set_atoms(atoms=get_atoms_from_file(filename))

        self.solvent = str(solvent) if solvent is not None else None        #: (str) Name of the solvent
        self.charge = int(charge)                                           #: (int) Charge in e
        self.mult = int(mult)                                               #: (int) Spin multiplicity 2S+1

        self.energy = None                                                  #: (float) Energy in Hartrees (Ha)

        self.reasonable_geometry = True                                     #: (bool)


class Molecule(BaseStruct):

    def get_charges(self, estimate=False):
        """
        Get the partial atomic charges using either XTB or estimate with RDKit using the Gasteiger charge scheme

        :param estimate: (bool)
        :param guess: (bool)
        :return:
        """

        if estimate and self.mol_obj is None:
            raise CgbindCritical('Cannot estimate charges without a rdkit molecule object')

        if estimate:
            try:
                rdPartialCharges.ComputeGasteigerCharges(self.mol_obj)
                charges = [float(self.mol_obj.GetAtomWithIdx(i).GetProp('_GasteigerCharge')) for i in range(self.n_atoms)]
            except:
                logger.error('RDKit failed to generate charges')
                return None

        else:
            charges = calculations.get_charges(self)

        return charges

    def _init_smiles(self, smiles, use_etdg_confs=False):
        """
        Initialise a Molecule object from a SMILES sting using RDKit
        :param smiles: (str) SMILES string
        :param use_etdg_confs: (bool) override the default conformer generation and use the ETDG algorithm
        :return:
        """
        logger.info('Initialising a Molecule from a SMILES string')
        try:
            self.mol_obj = Chem.MolFromSmiles(smiles)
            self.mol_obj = Chem.AddHs(self.mol_obj)
            self.charge = Chem.GetFormalCharge(self.mol_obj)
            self.n_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(self.mol_obj)
            self.n_h_donors = rdMolDescriptors.CalcNumHBD(self.mol_obj)
            self.n_h_acceptors = rdMolDescriptors.CalcNumHBA(self.mol_obj)

        except:
            logger.error('RDKit failed to generate mol objects')
            return

        logger.info('Running conformation generation with RDKit... running')
        method = AllChem.ETKDGv2() if use_etdg_confs is False else AllChem.ETDG()
        method.pruneRmsThresh = 0.3
        method.numThreads = Config.n_cores
        conf_ids = list(AllChem.EmbedMultipleConfs(self.mol_obj, numConfs=self.n_confs, params=method))
        logger.info('                                          ... done')

        try:
            self.volume = AllChem.ComputeMolVolume(self.mol_obj)
        except ValueError:
            logger.error('RDKit failed to compute the molecular volume')
            return

        self.bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in self.mol_obj.GetBonds()]
        self.conformers = extract_conformers_from_rdkit_mol_object(mol_obj=self.mol_obj, conf_ids=conf_ids)

        # Default to the first generated conformer in the absence of any other information
        self.set_atoms(atoms=self.conformers[0].atoms)

        return None

    def __init__(self, smiles=None, name='molecule', charge=0, mult=1, n_confs=1, filename=None, solvent=None,
                 use_etdg_confs=False):
        """
        Molecule. Inherits from cgbind.molecule.BaseStruct

        :param smiles: (str) SMILES string
        :param name: (str) Molecule name
        :param n_confs: (int) Number of conformers to initialise with
        :param charge: (int) Charge on the molecule
        :param mult: (int) Spin multiplicity on the molecule
        :param filename: (str)
        :param use_etdg_confs: (bool) Use an alternate conformer generation algorithm
        """
        logger.info('Initialising a Molecule object for {}'.format(name))

        super(Molecule, self).__init__(name=name, charge=charge, mult=mult, filename=filename, solvent=solvent)

        self.smiles = smiles                                        #: (str) SMILES string
        self.n_confs = n_confs                                      #: (int) Number of conformers initialised with

        self.mol_obj = None                                         #: (RDKit.mol object)

        self.n_rot_bonds = None                                     #: (int) Number of rotatable bonds
        self.n_h_donors = None                                      #: (int) Number of H-bond donors
        self.n_h_acceptors = None                                   #: (int) Number of H-bond acceptors
        self.volume = None                                          #: (float) Molecular volume in Å^3
        self.bonds = None                                           #: (list(tuple)) List of bonds defined by atom ids

        self.conformers = None                                      #: (list(BaseStruct)) List of conformers

        if smiles:
            self._init_smiles(smiles, use_etdg_confs=use_etdg_confs)

        if filename is not None:
            self.bonds = get_bond_list_from_atoms(self.atoms)
            self.conformers = [deepcopy(self)]
