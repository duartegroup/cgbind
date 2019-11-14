from cgbind.log import logger
from cgbind.molecule import Molecule


class Substrate(Molecule):

    def __init__(self, smiles=None, name='substrate', n_confs=10, charge=0, xyzs=None, conf_id=0):
        logger.info('Initialising a Substrate object for {}'.format(name))
        super(Substrate, self).__init__(smiles=smiles, name=name, charge=charge, n_confs=n_confs, xyzs=xyzs)

        logger.warning(f'Choosing the {conf_id} conformer of {self.name}')
        if self.conf_xyzs is not None:
            self.xyzs = self.conf_xyzs[conf_id]
