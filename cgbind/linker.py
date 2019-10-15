from cgbind.log import logger
from cgbind import m2l4
from cgbind import m4l6
from cgbind.molecule import Molecule


class Linker(Molecule):

    def set_xyzs_best_linker_conformer(self, conf_id):
        if conf_id is not None:
            self.xyzs = self.conf_xyzs[conf_id]
        else:
            logger.error("Didn't find a suitable conformer for {}".format(self.name))

    def __init__(self, smiles=None, name='linker', charge=0, n_confs=200, xyzs=None):

        logger.info('Initialising a Linker object for {}'.format(name))
        super(Linker, self).__init__(smiles=smiles, name=name, charge=charge, n_confs=n_confs, xyzs=xyzs)

        self.arch = None
        self.x_atom_ids = None   # The EXE for M2L4 and the XEEX atoms for M4L6

        # if self.arch == M2L4:
        #     self.exe_motifs = None   # Set in get_best_conformer
        #     self.conf_id = m2l4.get_best_conformer(self)

        #     if self.exe_motifs is not None:
        #         self.x_atom_ids = [i for motif in self.exe_motifs for i in motif]
        #         self.set_xyzs_best_linker_conformer(conf_id=self.conf_id)

        # if self.arch == M4L6:
        #     self.xeex_motifs = None   # Set in get_best_conformer
        #     self.conf_id = m4l6.get_best_linker_conformer(self)

        #     if self.xeex_motifs is not None:
        #         self.x_atom_ids = [i for motif in self.xeex_motifs for i in motif]
        #         self.set_xyzs_best_linker_conformer(conf_id=self.conf_id)
