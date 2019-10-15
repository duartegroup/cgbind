from cgbind.log import logger
import numpy as np
import itertools
from copy import deepcopy
from cgbind.molecule import Molecule
from cgbind.atoms import heteroatoms
from cgbind.geom import xyz2coord
from cgbind.templates import get_template
from cgbind.x_motifs import find_x_motifs
from cgbind.x_motifs import check_x_motifs
from cgbind.build import get_template_fitted_coords_and_cost


def get_shifted_template_x_motif_coords(linker_template, dr):
    """
    For a linker template modify the x motif coordinates by a particular distance (dr) along the shift vector

    e.g. for M2L4

               ^
    M -- X--   | shift vec
           |
           |
           |
    M -- X--   | shift vec


    :param linker_template:
    :param dr: (float) Distance in Å to shift the x motifs by
    :return:
    """

    shifted_coords = []

    for motif in linker_template.x_motifs:
        for coord in deepcopy(motif.coords):
            coord += dr * motif.norm_shift_vec
            shifted_coords.append(coord)

    return shifted_coords


class Linker(Molecule):

    def set_xyzs_best_linker_conformer(self, conf_id):
        if conf_id is not None:
            self.xyzs = self.conf_xyzs[conf_id]
        else:
            logger.error("Didn't find a suitable conformer for {}".format(self.name))

    def _find_possible_donor_atoms(self):
        return [i for i in range(self.n_atoms) if self.xyzs[i][0] in heteroatoms]

    def _strip_possible_x_motifs_on_connectivity(self):
        """
        For a list of x motifs remove those which don't have the same number of atoms as the template linkers
        :return: (list) new list of x motifs
        """
        return [x_motif for x_motif in self.x_motifs
                if x_motif.n_atoms == self.cage_template.linkers[0].x_motifs[0].n_atoms]

    def __init__(self, smiles=None, name='linker', charge=0, n_confs=200, xyzs=None, arch='m2l4'):

        logger.info('Initialising a Linker object for {}'.format(name))
        super(Linker, self).__init__(smiles=smiles, name=name, charge=charge, n_confs=n_confs, xyzs=xyzs)

        self.arch = arch
        self.cage_template = get_template(arch_name=arch)

        self.coords = xyz2coord(self.xyzs)
        self.centroid = np.average(self.coords, axis=0)

        self.x_atoms = self._find_possible_donor_atoms()
        self.x_motifs = find_x_motifs(self)
        check_x_motifs(self)
        self.x_motifs = self._strip_possible_x_motifs_on_connectivity()

        # For all the possible combinations of x_motifs minimise the RMSD between the x_motifs and the template
        # x_motifs. The template needs to be modified to accommodate longer linkers with the same architecture

        template_linker = self.cage_template.linkers[0]
        n_x_motifs_in_linker = len(template_linker.x_motifs)

        for x_motifs in itertools.combinations(self.x_motifs, n_x_motifs_in_linker):

            x_coords = [self.coords[atom_id] for motif in x_motifs for atom_id in motif.atom_ids]
            shifted_coords = get_shifted_template_x_motif_coords(linker_template=template_linker, dr=0.1)
            _, cost = get_template_fitted_coords_and_cost(self, template_x_coords=shifted_coords,
                                                          coords_to_fit=x_coords)

            print(cost)



            pass

        self.x_motifs = None
