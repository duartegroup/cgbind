from cgbind.log import logger
from cgbind.molecule import Molecule
from cgbind.atoms import heteroatoms
from rdkit.Chem import AllChem
from cgbind.confomers import extract_xyzs_from_rdkit_mol_object


class Substrate(Molecule):

    def gen_confs(self, n_confs=1):
        """Populate self.conf_xyzs by calling RDKit"""

        if self.smiles is None:
            logger.error('Could not generate conformers. Substrate was not initialised from a SMILES string')
            return None

        if self.mol_obj is None:
            logger.error('Could not generate conformers. Molecule did not have an associated RDKit mol_obj')
            return None

        self.conf_ids = list(AllChem.EmbedMultipleConfs(self.mol_obj, numConfs=n_confs, params=AllChem.ETKDG()))
        self.conf_xyzs = extract_xyzs_from_rdkit_mol_object(mol_obj=self.mol_obj, conf_ids=self.conf_ids)

        return None

    def __init__(self, smiles=None, name='substrate', n_confs=1, charge=0, mult=1, xyzs=None):
        logger.info('Initialising a Substrate object for {}'.format(name))
        super(Substrate, self).__init__(smiles=smiles, name=name, charge=charge, n_confs=n_confs, mult=mult, xyzs=xyzs)

        if self.n_atoms is None:
            logger.error('Molecule had no atoms')
            return

        self.x_atom_ids = [i for i in range(self.n_atoms) if self.xyzs[i][0] in heteroatoms]
