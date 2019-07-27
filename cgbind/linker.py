from .log import logger
from . import m2l4
from . import m4l6
from .architectures import M2L4
from .architectures import M4L6
from .architectures import archs
from .molecule import Molecule


class Linker(Molecule):

    def set_arch(self, arch):
        if arch not in archs:
            logger.critical('Could not set the architecture')
            exit()
        else:
            self.arch = arch

    def set_xyzs_best_linker_conformer(self, conf_id):
        if conf_id is not None:
            self.xyzs = self.conf_xyzs[conf_id]
        else:
            logger.error("Didn't find a suitable conformer for {}".format(self.name))

    def __init__(self, smiles=None, name='linker', charge=0, n_confs=200, arch=M2L4, xyzs=None):

        logger.info('Initialising a Linker object for {}'.format(name))
        super(Linker, self).__init__(smiles=smiles, name=name, charge=charge, n_confs=n_confs, xyzs=xyzs)

        self.arch = None
        self.set_arch(arch)

        if self.arch == M2L4:
            self.x_atom_ids = None    # Set in get_best_linker_conformer
            self.x_x_midpoint = None  # Set in get_mid_atom

            self.conf_id = m2l4.get_best_linker_conformer(self)
            self.set_xyzs_best_linker_conformer(conf_id=self.conf_id)
            self.mid_atom_id = m2l4.get_mid_atom_id(self)

        if self.arch == M4L6:
            self.xeex_motifs = None   # Set in get_best_conformer
            self.conf_id = m4l6.get_best_linker_conformer(self)
            self.set_xyzs_best_linker_conformer(conf_id=self.conf_id)
