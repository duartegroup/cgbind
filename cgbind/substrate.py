from cgbind.log import logger
from cgbind.molecule import Molecule
from rdkit.Chem import AllChem
from cgbind.molecule import extract_conformers_from_rdkit_mol_object


class Substrate(Molecule):

    def gen_confs(self, n_confs=1):
        """Populate self.conf_xyzs by calling RDKit"""

        if self.smiles is None:
            logger.error('Could not generate conformers. Substrate was not initialised from a SMILES string')
            return None

        if self.mol_obj is None:
            logger.error('Could not generate conformers. Molecule did not have an associated RDKit mol_obj')
            return None

        conf_ids = list(AllChem.EmbedMultipleConfs(self.mol_obj, numConfs=n_confs, params=AllChem.ETKDG()))
        self.conformers = extract_conformers_from_rdkit_mol_object(mol_obj=self.mol_obj, conf_ids=conf_ids)

        return None

    def __init__(self, smiles=None, name='substrate', n_confs=1, charge=0, mult=1, filename=None, solvent=None):
        """
        Substrate. Inherits from cgbind.molecule.Molecule

        :param smiles: (str) SMILES string
        :param name: (str) Molecule name
        :param n_confs: (int) Number of conformers to initialise with
        :param charge: (int) Charge on the molecule
        :param mult: (int) Spin multiplicity on the molecule
        :param filename: (str)
        """

        logger.info('Initialising a Substrate object for {}'.format(name))
        super(Substrate, self).__init__(smiles=smiles, name=name, charge=charge, n_confs=n_confs,
                                        mult=mult, filename=filename, solvent=solvent)
