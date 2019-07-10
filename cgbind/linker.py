import numpy as np
from .log import logger
from . import m2l4
from . import m4l6
from .config import Config
from .optimisation import opt_geom
from rdkit.Chem import AllChem
from rdkit import Chem
from .confomers import gen_conformer_mol_files
from .confomers import confomer_mol_files_to_xyzs
from .geom import xyz2coord
from .geom import calc_midpoint


class Linker(object):

    def _get_mid_atom_id(self):
        """
        For a linker then the 'mid_atom' is needed to define a plane from which the M2L4 cage may be built. This is
        found by calculating the abs difference in the A-X1, A-X2 distances for all atoms A, the minimum of which will
        be the atom closest to the midpoint of the linker as possible.
        :return: The ID of the atom in the xyz list
        """
        logger.info('Calculating mid atom id for linker {}'.format(self.name))

        mid_atom_id = 0
        x1_atom_id, x2_atom_id = self.x_atom_ids
        x1_coords, x2_coords = xyz2coord(self.xyzs[x1_atom_id]), xyz2coord(self.xyzs[x2_atom_id])
        self.x_x_midpoint = calc_midpoint(x1_coords, x2_coords)

        min_dist_diff = 999.9

        for i in range(len(self.xyzs)):
            atom_i_coords = xyz2coord(self.xyzs[i])
            dist_diff = np.abs(np.linalg.norm(atom_i_coords - x1_coords) - np.linalg.norm(atom_i_coords - x2_coords))
            if dist_diff < min_dist_diff:
                min_dist_diff = dist_diff
                mid_atom_id = i

        return mid_atom_id

    def __init__(self, smiles=None, name='linker', charge=0, n_confs=200, opt=True, arch='m2l4', n_cores=1):
        """
        Initialise a Linker object
        :param smiles: SMILES string of the linker (str)
        :param name: Name of the linker (str)
        :param charge: Total charge on a single linker (int)
        :param n_confs: Number of conformers to generate. A too low number here may result in get_best_linker_conformer
        to not find the best conformer (int)
        :param opt: Optimise the linker geometry with the method specified in Config.opt_method
        :param arch: Cage architecture (str)
        :param n_cores: Number of cores to use to perform the geometry optimisation
        """
        logger.info('Initialising a Linker object for {}'.format(name))

        self.name = name
        self.smiles = smiles
        self.charge = charge
        self.arch = arch.lower()
        self.energy, self.xyzs = None, None

        if not smiles:
            logger.warning('Linker {} was not initialised with a SMILES string'.format(name))
            return

        try:
            self.mol_obj = Chem.MolFromSmiles(smiles)
            self.mol_obj = Chem.AddHs(self.mol_obj)
            self.charge = Chem.GetFormalCharge(self.mol_obj)

        except RuntimeError:
            logger.error('RDKit failed to generate mol objects')
            return

        if not Config.suppress_print:
            print("{:<30s}{:<50s}{:>10s}".format('Confomer generation for', self.name, 'Running'))
        self.conf_ids = list(AllChem.EmbedMultipleConfs(self.mol_obj, numConfs=n_confs, params=AllChem.ETKDG()))
        self.conf_filenames = [self.name + '_conf' + str(i) + '.mol' for i in self.conf_ids]
        gen_conformer_mol_files(self)

        self.n_atoms = self.mol_obj.GetNumAtoms()
        self.conf_xyzs = confomer_mol_files_to_xyzs(self.conf_filenames, self.n_atoms)
        self.x_atom_ids = []                                 # Set in get_best_conformer

        if self.arch == 'm2l4':
            self.conf_id = m2l4.get_best_linker_conformer(self)
            if self.conf_id is not None:
                self.xyzs = self.conf_xyzs[self.conf_id]
            else:
                logger.error("Didn't find a suitable conformer for {}".format(self.name))

            if self.xyzs:
                self.x_x_midpoint = None  # Set in get_mid_atom
                self.mid_atom_id = self._get_mid_atom_id()

        if self.arch == 'm4l6':
            self.xeex_motifs = None      # Set in get_best_conformer
            self.conf_id = m4l6.get_best_linker_conformer(self)
            if self.conf_id is not None:
                self.xyzs = self.conf_xyzs[self.conf_id]

        if opt:
            self.xyzs, self.energy = opt_geom(self.xyzs, self.name, charge=charge, n_cores=n_cores)
