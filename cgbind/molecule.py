from rdkit.Chem import AllChem
from rdkit import Chem
from cgbind.log import logger
from cgbind.optimisation import opt_geom
from cgbind.single_point import singlepoint
from cgbind.input_output import print_output
from cgbind.input_output import xyzs2xyzfile
from cgbind.confomers import gen_conformer_mol_files
from cgbind.confomers import confomer_mol_files_to_xyzs
from cgbind.geom import calc_com


class Molecule(object):

    def optimise(self, n_cores):
        self.xyzs, self.energy = opt_geom(self.xyzs, self.name, charge=self.charge, n_cores=n_cores)

    def singlepoint(self, n_cores=1):
        self.energy = singlepoint(self, n_cores)

    def print_xyzfile(self):
        xyzs2xyzfile(xyzs=self.xyzs, basename=self.name)

    def set_charge(self, charge):
        assert type(charge) == int
        self.charge = charge

    def set_com(self):
        self.com = calc_com(self.xyzs)

    def init_smiles(self, smiles):
        """
        Initialise a Molecule object from a SMILES sting using RDKit
        :param smiles: (str) SMILES string
        :return:
        """
        logger.info('Initialising a Molecule from a SMILES strings ')
        try:
            self.mol_obj = Chem.MolFromSmiles(smiles)
            self.mol_obj = Chem.AddHs(self.mol_obj)
            self.charge = Chem.GetFormalCharge(self.mol_obj)

        except RuntimeError:
            logger.error('RDKit failed to generate mol objects')
            return

        print_output('Confomer generation for', self.name, 'Running')
        self.conf_ids = list(AllChem.EmbedMultipleConfs(self.mol_obj, numConfs=self.n_confs, params=AllChem.ETKDG()))
        self.conf_filenames = [self.name + '_conf' + str(i) + '.mol' for i in self.conf_ids]
        gen_conformer_mol_files(self)
        print_output('', '', 'Done')

        self.n_atoms = self.mol_obj.GetNumAtoms()
        self.conf_xyzs = confomer_mol_files_to_xyzs(self.conf_filenames, self.n_atoms)
        self.xyzs = self.conf_xyzs[0]

    def __init__(self, smiles=None, name='molecule', charge=0, n_confs=10, xyzs=None):
        logger.info('Initialising a Molecule object for {}'.format(name))

        self.name = name
        self.smiles = smiles
        self.xyzs = xyzs
        self.n_confs = n_confs

        self.charge = None
        self.set_charge(charge)

        self.energy = None
        self.mol_obj = None
        self.n_atoms = None
        self.com = None

        self.conf_ids = None
        self.conf_filenames = None
        self.conf_xyzs = None

        if smiles:
            self.init_smiles(smiles)

        if xyzs:
            self.n_atoms = len(xyzs)
