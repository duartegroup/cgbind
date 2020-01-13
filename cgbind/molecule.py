from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdPartialCharges
from cgbind.log import logger
from cgbind.config import Config
from cgbind.input_output import print_output
from cgbind.input_output import xyzs2xyzfile
from cgbind.confomers import extract_xyzs_from_rdkit_mol_object
from cgbind.geom import calc_com
from autode.bond_lengths import get_bond_list_from_rdkit_bonds
from autode.bond_lengths import get_xyz_bond_list
from cgbind import calculations


class BaseStruct:

    def print_xyzfile(self, force=False):
        """
        Print a .xyz file from self.xyzs provided self.reasonable_geometry is True

        :param force: (bool) Force the printing of the .xyz if self.reasonable_geometry is False
        :return: None
        """

        if self.reasonable_geometry or force:
            xyzs2xyzfile(xyzs=self.xyzs, basename=self.name)

        return None

    def singlepoint(self, method, keywords=None, n_cores=1, max_core_mb=1000):
        """
        Perform a single-point energy evaluation using an electronic structure theory method e.g. XTB, ORCA, G09

        :param method: (autode.ElectronicStructureMethod)
        :param keywords: (list(str)) Keywords to use for the ESM e.g. ['SP', 'PBE', 'def2-SVP']
        :param n_cores: (int) Number of cores for the calculation to use
        :param max_core_mb: (float) Number of megabytes of memory per core to use e.g. n_cores=2: max_core_mb=4000 =>
                                    max 8 GB total memory usage
        :return: None
        """
        return calculations.singlepoint(self, method, keywords, n_cores, max_core_mb)

    def optimise(self, method, keywords=None, n_cores=1, max_core_mb=1000, cartesian_constraints=None):
        """
        Perform a single-point energy evaluation using an electronic structure theory method e.g. XTB, ORCA, G09

        :param method: (autode.ElectronicStructureMethod)
        :param keywords: (list(str)) Keywords to use for the ESM e.g. ['Opt', 'PBE', 'def2-SVP']
        :param n_cores: (int) Number of cores for the calculation to use
        :param max_core_mb: (float) Number of megabytes of memory per core to use
        :param cartesian_constraints: (list(int)) List of atom ids to constrain to their current coordinates
        :return: None
        """
        return calculations.optimise(self, method, keywords, n_cores, max_core_mb, cartesian_constraints)

    def set_xyzs(self, xyzs):
        """
        Set the xyzs of a molecular structure

        :param xyzs: (list(list)) e.g [[C, 0.0, 0.0, 0.0], ...]
        :return: None
        """

        if xyzs is not None:
            assert type(xyzs) == list
            assert type(xyzs[0]) == list
            assert len(xyzs[0]) == 4

            self.xyzs = xyzs
            self.n_atoms = len(xyzs)

            logger.info('Successfully set xyzs and n_atoms')

        else:
            self.reasonable_geometry = False
            logger.warning('xyzs were None -> n_atoms also None & geometry is *not* reasonable')

        return None

    def __init__(self, name='molecule', charge=0, mult=1, xyzs=None, solvent=None):
        """
        Base structure class

        :param name: (str)
        :param charge: (int)
        :param mult: (int) Spin multiplicity
        :param xyzs: (list(list))
        :param solvent: (str)
        """

        self.name = str(name)                                               #: (str) Name of the structure

        self.n_atoms = None                                                 #: (int) Number of atoms
        self.xyzs = None                                                    #: (list(list)) Geometry of the structure
        self.set_xyzs(xyzs)

        self.solvent = str(solvent) if solvent is not None else None        #: (str) Name of the solvent
        self.charge = int(charge)                                           #: (int) Charge in e
        self.mult = int(mult)                                               #: (int) Spin multiplicity 2S+1

        self.energy = None                                                  #: (float) Energy in Hartrees (Ha)

        self.reasonable_geometry = True                                     #: (bool)


class Molecule(BaseStruct):

    def set_com(self):
        self.com = calc_com(self.xyzs)

    def get_charges(self, estimate=False, guess=False):
        """
        Get the partial atomic charges using either XTB or estimate with RDKit using the Gasteiger charge scheme

        :param estimate: (bool)
        :param guess: (bool)
        :return:
        """

        charges = None

        if estimate:
            rdPartialCharges.ComputeGasteigerCharges(self.mol_obj)
            try:
                charges = [float(self.mol_obj.GetAtomWithIdx(i).GetProp('_GasteigerCharge')) for i in range(self.n_atoms)]
            except:
                logger.error('RDKit failed to generate charges')
                return None

        elif guess:
            # TODO write this function
            pass

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
        logger.info('Initialising a Molecule from a SMILES strings ')
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

        print_output('Conformer generation for', self.name, 'Running')
        method = AllChem.ETKDG() if use_etdg_confs is False else AllChem.ETDG()
        method.pruneRmsThresh = 0.3
        method.numThreads = Config.n_cores
        conf_ids = list(AllChem.EmbedMultipleConfs(self.mol_obj, numConfs=self.n_confs, params=method))
        try:
            self.volume = AllChem.ComputeMolVolume(self.mol_obj)
        except ValueError:
            logger.error('RDKit failed to compute the molecular volume')
            return

        self.bonds = get_bond_list_from_rdkit_bonds(rdkit_bonds_obj=self.mol_obj.GetBonds())
        print_output('', '', 'Done')

        self.conf_xyzs = extract_xyzs_from_rdkit_mol_object(mol_obj=self.mol_obj, conf_ids=conf_ids)
        self.set_xyzs(xyzs=self.conf_xyzs[0])

        return None

    def __init__(self, smiles=None, name='molecule', charge=0, mult=1, n_confs=1, xyzs=None, solvent=None,
                 use_etdg_confs=False):
        """
        Molecule. Inherits from cgbind.molecule.BaseStruct

        :param smiles: (str) SMILES string
        :param name: (str) Molecule name
        :param n_confs: (int) Number of conformers to initialise with
        :param charge: (int) Charge on the molecule
        :param mult: (int) Spin multiplicity on the molecule
        :param xyzs: (list(list))
        :param use_etdg_confs: (bool) Use an alternate conformer generation algorithm
        """
        logger.info('Initialising a Molecule object for {}'.format(name))

        super(Molecule, self).__init__(name=name, charge=charge, mult=mult, xyzs=xyzs, solvent=solvent)

        self.smiles = smiles                                        #: (str) SMILES string
        self.n_confs = n_confs                                      #: (int) Number of conformers initialised with

        self.mol_obj = None                                         #: (RDKit.mol object)
        self.com = None                                             #: (np.ndarray) Center of mass (x, y, z)

        self.n_rot_bonds = None                                     #: (int) Number of rotatable bonds
        self.n_h_donors = None                                      #: (int) Number of H-bond donors
        self.n_h_acceptors = None                                   #: (int) Number of H-bond acceptors
        self.volume = None                                          #: (float) Molecular volume in Å^3
        self.bonds = None                                           #: (list(tuple)) List of bonds defined by atom ids

        self.conf_xyzs = None                                       #: (list(xyzs)) List of xyzs of the conformers

        if smiles:
            self._init_smiles(smiles, use_etdg_confs=use_etdg_confs)

        if xyzs is not None:
            self.bonds = get_xyz_bond_list(xyzs=self.xyzs)
        else:
            logger.error('Failed to generate or set molecular xyzs. self.bonds is None')
