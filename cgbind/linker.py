from cgbind.log import logger
import numpy as np
import itertools
from cgbind.molecule import Molecule
from cgbind.atoms import heteroatoms
from cgbind.geom import xyz2coord
from cgbind.templates import get_template
from cgbind.x_motifs import find_x_motifs
from cgbind.x_motifs import check_x_motifs
from cgbind.build import get_template_fitted_coords_and_cost


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
        return [x_motif for x_motif in self.x_motifs if len(x_motif) == self.cage_template.linkers[0].len_x_motif]

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
        min_cost, best_x_motifs = None, None

        def get_x_motifs_coords_new_dist(motifs_coords, motifs_vecs, r, curr_rs):

            for i, coords in enumerate(motifs_coords):
                coords += (r - curr_rs[i]) * motifs_vecs[i] / curr_rs[i]

            return [coord for coords in motifs_coords for coord in coords]

        for x_motifs in itertools.combinations(self.x_motifs, n_x_motifs_in_linker):

            template_x_motifs_coords = [np.array([template_linker.coords[i] for i in motif])
                                        for motif in template_linker.x_motifs]

            shift_vecs = [np.average(coords, axis=0) - template_linker.centroid for coords in template_x_motifs_coords]
            curr_rs = [np.linalg.norm(vec) for vec in shift_vecs]

            template_coords = get_x_motifs_coords_new_dist(motifs_coords=template_x_motifs_coords,
                                                           motifs_vecs=shift_vecs, r=7.0, curr_rs=curr_rs)
            coords = [self.coords[i] for motif in x_motifs for i in motif]

            _, cost = get_template_fitted_coords_and_cost(self, template_coords, coords_to_fit=coords)

            print(cost)
            exit()




            pass

        self.x_motifs = best_x_motifs


