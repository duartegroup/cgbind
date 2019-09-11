import numpy as np
from cgbind.log import logger
from cgbind.config import Config
from cgbind.atoms import heteroatoms
from cgbind.molecule import Molecule


def get_bond_vector_angles_and_x_ids(xyzs):

        ideal_angle = np.pi                         # Ideal angle between X- and -X vectors for binding

        # ---------------------------- Get x(hetero)atom ids -----------------------
        x_ids = []
        for i in range(len(xyzs)):
            if xyzs[i][0] in heteroatoms:
                x_ids.append(i)

        # ------------------------ Get x(hetero)atom bond vectors -----------------------
        x_bond_vectors = []
        for i in x_ids:
            min_dist = 999.9
            bond_vec = np.zeros(3)
            for j in range(len(xyzs)):
                vec = np.array(xyzs[i][1:4]) - np.array(xyzs[j][1:4])
                dist = np.linalg.norm(vec)
                if dist < min_dist and i != j:
                    bond_vec = vec
                    min_dist = dist

            x_bond_vectors.append(bond_vec)

        # ----------- For a suitable substrate the vectors need to be close to inverses ------
        closest_angle = 0
        closest_x_ids = [None, None]

        for i in range(len(x_bond_vectors)):
            for j in range(len(x_bond_vectors)):
                if i > j:
                    norm_vec_i = x_bond_vectors[i] / np.linalg.norm(x_bond_vectors[i])
                    norm_vec_j = x_bond_vectors[j] / np.linalg.norm(x_bond_vectors[j])
                    angle_rad = np.arccos(np.dot(norm_vec_i, norm_vec_j))

                    if np.abs(angle_rad - ideal_angle) < np.abs(closest_angle - ideal_angle):

                        # --------- Distances need to be reasonable i.e. 5Å < dist < 6Å ------------

                        dist = np.linalg.norm(np.array(xyzs[x_ids[i]][1:4]) - np.array(xyzs[x_ids[j]][1:4]))

                        if 4 < dist < 8:
                            closest_angle = angle_rad
                            closest_x_ids = [x_ids[i], x_ids[j]]

        return closest_angle, closest_x_ids


class Substrate(Molecule):

    def get_conformer_with_collinear_x_x(self):
        """
        This function is set up to return the 'best' confomer of a substrate suitable for binding in an M2L4 cage
        no consideration of the energy is made. The user should check that the confomer that is output of this function
        is accessible with kT. Furthermore, many conformers may contribute to the binding affinity; this is also not
        considered (only one conformer is selected).

        :return: confomer ids
        """
        logger.info('Getting conformer with collinear_X_X for {}'.format(self.name))

        best_conf_id = None
        best_x_ids = [None, None]
        best_angle = 0
        ideal_angle = np.pi                         # Ideal angle between X- and -X vectors for binding

        for conf_id in self.conf_ids:
            closest_angle, closest_x_ids = get_bond_vector_angles_and_x_ids(xyzs=self.conf_xyzs[conf_id])

            if np.abs(closest_angle - ideal_angle) < np.abs(best_angle - ideal_angle):
                best_angle = closest_angle
                best_conf_id = conf_id
                best_x_ids = closest_x_ids

        self.x_x_atom_ids = best_x_ids

        if not best_conf_id and not Config.suppress_print:
            print("Didn't find a conformer with collinear E-X")

        return best_conf_id

    def get_n_heteroatoms(self):
        """
        Get the number of heteroatoms (X) in a substrate. If >1 then will check if we can align two X-E with the M-M
        vector, if =1 then will try and align the X-E bond vector with the M-M vector
        :return: Number of heteroatoms
        """
        logger.info('Getting the number of heteroatoms for {}'.format(self.name))
        if self.xyzs is None:
            return 0

        n_hetero_atoms = 0
        for xyz_line in self.xyzs:
            if xyz_line[0] in heteroatoms:
                n_hetero_atoms += 1
        return n_hetero_atoms

    def get_x_atom_ids(self):
        logger.info('Getting heteroatom ids in substrate')
        return [i for i in range(len(self.xyzs)) if self.xyzs[i][0] in heteroatoms]

    def __init__(self, smiles=None, name='substrate', n_confs=10, charge=0, xyzs=None):
        logger.info('Initialising a Substrate object for {}'.format(name))
        super(Substrate, self).__init__(smiles=smiles, name=name, charge=charge, n_confs=n_confs, xyzs=xyzs)
        self.x_x_atom_ids = None
        self.x_atom_ids = None

        self.n_heteroatoms = self.get_n_heteroatoms()
        logger.info('Substrate has {} heteroatoms'.format(self.n_heteroatoms))

        if self.n_heteroatoms >= 2:
            self.conf_id_collinear_XX = self.get_conformer_with_collinear_x_x()
            if self.conf_id_collinear_XX and self.conf_xyzs is not None:
                # If there is a conformer with collinear heteroatoms the there is just one binding mode in M2L4
                self.xyzs = self.conf_xyzs[self.conf_id_collinear_XX]

        elif self.n_heteroatoms > 0:
            # Choose the first conformer. This may not the lowest energy, but the lowest E for the free ligand might
            # not be the one bound inside the cage, for the optimum binding affinity -- this is a hard problem
            logger.warning('Choosing the first conformer of {}'.format(self.name))
            if self.conf_xyzs is not None:
                self.xyzs = self.conf_xyzs[0]
            self.x_atom_ids = self.get_x_atom_ids()

        else:
            # Choose the first conformer, again perhaps not a great idea
            logger.warning('Choosing the first conformer of {}'.format(self.name))
            if self.conf_xyzs is not None:
                self.xyzs = self.conf_xyzs[0]
