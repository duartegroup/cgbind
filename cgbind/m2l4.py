import numpy as np
import itertools
from copy import deepcopy
from cgbind.log import logger
from cgbind.geom import rotation_matrix
from cgbind.geom import xyz2coord
from cgbind.geom import get_closest_bonded_atom_id
from cgbind.geom import calc_normalised_vector
from cgbind.geom import calc_com
from cgbind.geom import calc_cdist
from cgbind.geom import calc_midpoint
from cgbind.geom import cat_cage_subst_coords
from cgbind.build import get_cost_xyzs
from scipy.optimize import minimize


def build(cage, linker):
    """
    From a linker object construct a M2L4 metallocage by fitting to the template as a function of the M-M distance
    to minimise the cost function


            N
           /|
    --N - M- N--
    |    /  |  |
    |   N   |  |
    |  |    |  |
    |  |    N  |
    |  |   /   |
    --N|- M- N--
       |/
       N


    :param cage: A Cage object
    :param linker: A linker object
    """
    logger.info('Building M2L4 metallocage'.format(cage.name))

    if any(attr is None for attr in (linker.xyzs, linker.x_atom_ids)):
        logger.error('A required linker property was None')
        return None

    opt_res = minimize(get_cost_xyzs, x0=np.array([5.0]),  method='BFGS',
                       args=(cage, linker, Template(), get_xx_coords(linker), False))

    cost, xyzs = get_cost_xyzs(r=opt_res.x[0], cage=cage, linker=linker, template=Template(),
                               coords_to_fit=get_xx_coords(linker), return_xyzs=True)

    if len(xyzs) != cage.arch.n_metals + cage.arch.n_linkers * linker.n_atoms:
        logger.error('Could not generate M2L4 metallocage')
    else:
        logger.info('Built M2L4 metallocage')
        cage.xyzs = xyzs

    return None


def get_xx_coords(linker):
    return [xyz2coord(xyzs=linker.xyzs[i]) for i in range(linker.n_atoms) if i in linker.x_atom_ids]


def get_best_conformer(linker, enclosed_cavity_rad=2.0):
    """
    For a linker to be suitable for use in a M2L4 cage we need co-linear N-M vectors, which are approx the N-'lone pair'
    vectors

    ------N -->
    |
    |
    |     empty
    |
    |
    ------N -->


    Also, there needs to be a vacant cavity between the two nitrogen atoms, as above.

    :param linker: A linker object
    :param enclosed_cavity_rad:
    :return: Best conformer id
    """
    logger.info('Getting the best linker conformer for {}'.format(linker.name))
    logger.warning('Will only consider N atoms as donors')

    n_atom_ids = get_donor_nitrogen_atom_ids(linker=linker)

    if linker.conf_ids is None or linker.conf_xyzs is None or n_atom_ids is None:
        logger.error('A required linker property was None')
        return None

    # For each pair of N atoms find the conformer with the most co-linear N-lp vector and has no atoms within
    # a sphere defined by enclosed_cavity_rad
    ideal_cos_theta = 1.0                                       # For co-linear vectors theta=0 and cos(0) = 1
    curr_cos_theta = 9999.9
    best_conf_id = None

    for (n1, n2) in itertools.combinations(n_atom_ids, 2):
        for conf_id in linker.conf_ids:

            lp_vectors = get_nitrogen_lone_pair_vectors(linker=linker, nitrogen_atom_ids=n_atom_ids, conf_id=conf_id)
            cos_theta = np.dot(lp_vectors[0], lp_vectors[1])

            # If the angle is smaller then we have a better conformer
            if np.abs(cos_theta - ideal_cos_theta) < np.abs(curr_cos_theta - ideal_cos_theta) or best_conf_id is None:
                curr_cos_theta = cos_theta
                best_conf_id = conf_id

        # Check that there is a vacant cavity suitable to construct a M2L4 cage
        linker_coords = xyz2coord(linker.conf_xyzs[best_conf_id])
        n_n_midpoint = calc_midpoint(coord1=linker_coords[n1], coord2=linker_coords[n2])
        atom_midpoint_dists = [np.linalg.norm(n_n_midpoint - linker_coords[i]) for i in range(linker.n_atoms)]

        linker_has_enclosed_cavity = all([dist > enclosed_cavity_rad for dist in atom_midpoint_dists])

        if linker_has_enclosed_cavity:
            logger.info('Found conformer with vacant cavity suitable to construct a M2L4 cage')
            populate_exe_motifs(linker=linker, n_id1=n1, n_id2=n2)

            logger.warning('If there are multiple possibilities for reasonable donor atoms they will be excluded')
            return best_conf_id

    return None


def populate_exe_motifs(linker, n_id1, n_id2):
    logger.info('Populating EXE motifs of the linker with the corresponding atom ids')
    linker.exe_motifs = []

    n_n_midpoint = calc_midpoint(coord1=xyz2coord(linker.xyzs[n_id1]), coord2=xyz2coord(linker.xyzs[n_id2]))
    for n_atom_id in [n_id1, n_id2]:
        atoms_bonded_to_n = [[atom_id for atom_id in bond if atom_id != n_atom_id][0]
                             for bond in linker.bonds if n_atom_id in bond]

        assert len(atoms_bonded_to_n) == 2
        e1, e2 = atoms_bonded_to_n

        # Ensure the ordering is correct to fit the two EXE motifs of each linker
        if calc_cdist(xyz2coord(linker.xyzs[e1]), n_n_midpoint) > calc_cdist(xyz2coord(linker.xyzs[e2]), n_n_midpoint):
            exe_motif = [e1, n_atom_id, e2]
        else:
            exe_motif = [e2, n_atom_id, e1]

        linker.exe_motifs.append(exe_motif)

    if len(linker.exe_motifs) == 0:
        logger.error('Could not find and EXE motifs')
        linker.exe_motifs = None

    return None


def get_nitrogen_lone_pair_vectors(linker, nitrogen_atom_ids, conf_id=None):
    """

    :param linker: linker object
    :param nitrogen_atom_ids: (list(int)) list of donor nitrogen ids
    :param conf_id: (int) id of the conformer to use in getting the xyzs
    :return: (dict) keyed with the nitrogen atom id and valued with the normalised N-lp vector (np.ndarray)
    """

    lp_vectors = []
    linker.exe_motifs = []
    assert len(nitrogen_atom_ids) == 2

    for nitrogen_atom_id in nitrogen_atom_ids:
        bonded_atom_ids = [[atom_id for atom_id in bond if atom_id != nitrogen_atom_id][0]
                           for bond in linker.bonds if nitrogen_atom_id in bond]

        if conf_id is not None:
            xyzs = linker.conf_xyzs[conf_id]
        else:
            xyzs = linker.xyzs

        # Get the vector that corresponds to the 'lone pair' of the N, which will be negative of the average bond vec
        lp_vec = np.sum(np.array([(xyz2coord(xyzs[nitrogen_atom_id]) - xyz2coord(xyzs[atom_id]))
                                  for atom_id in bonded_atom_ids]), axis=0) / float(len(bonded_atom_ids))

        lp_vectors.append(lp_vec / np.linalg.norm(lp_vec))

    return lp_vectors


def get_donor_nitrogen_atom_ids(linker):
    """
    :param linker: linker object
    :return: (list) nitrogen atom ids
    """

    if linker.n_atoms is None or linker.xyzs is None:
        logger.error('A required linker property was None')
        return None

    nitrogen_atom_ids = [i for i in range(linker.n_atoms) if linker.xyzs[i][0] == 'N']
    if len(nitrogen_atom_ids) < 2:
        logger.error('There are fewer than 2 nitrogen atoms in the linker. Cannot build a cage')
        return None

    for n_atom_id in nitrogen_atom_ids.copy():
        n_bonds_to_n = len([bond for bond in linker.bonds if n_atom_id in bond])

        # If the number of bonds to the nitrogen atom is 2 we have a 'lone-pair' and can form a cage
        if n_bonds_to_n != 2:
            nitrogen_atom_ids.remove(n_atom_id)

    return nitrogen_atom_ids


def rot_minimise_repulsion(subst_coords, cage_coords, rot_axes, n_rot_steps=100):
    """
    Rotate a substrate around the M-M (z) axis as to reduce the steric repulsion between it and the cage
    :param subst_coords:
    :param cage_coords:
    :param rot_axes: List of orthogonal normalised vectors
    :param n_rot_steps: Number of rotation steps to perform in each axis
    :return:
    """
    logger.info('Minimising steric repulsion between substrate and cage by rotation')

    best_theta = 0.0
    best_sum_cage_substrate_inverse_dists = 9999999.9
    best_axis = rot_axes[0]                         # TODO allow for a combination of axes

    for rot_axis in rot_axes:
        for theta in np.linspace(0, np.pi / 2.0, n_rot_steps):
            tmp_substrate_coords = deepcopy(subst_coords)
            rot_matrix = rotation_matrix(rot_axis, theta)
            tmp_rot_substrate_coords = [np.matmul(rot_matrix, coord) for coord in tmp_substrate_coords]

            sum_cage_substrate_inverse_dists = 0.0
            for cage_coord in cage_coords:
                for substrate_coord in tmp_rot_substrate_coords:
                    dist = np.linalg.norm(cage_coord - substrate_coord)
                    if dist < 2.0:
                        sum_cage_substrate_inverse_dists += 100.0 / dist  # Penalise very short distances

            if sum_cage_substrate_inverse_dists < best_sum_cage_substrate_inverse_dists:
                best_theta = theta
                best_sum_cage_substrate_inverse_dists = sum_cage_substrate_inverse_dists
                best_axis = rot_axis

    rot_matrix = rotation_matrix(best_axis, best_theta)
    rot_substrate_coords = [np.matmul(rot_matrix, coord) for coord in subst_coords]

    return rot_substrate_coords


def add_substrate_x_x(cage, substrate):
    """
    Add the substrate xyzs to the list of cage xyzs. Locate so the X-X vector is aligned with the M-M vector,
    then rotate in the M-M axis to minimise steric clash between the substrate and the cage

              N
             /|
    ----N - M- N----
    |      /  |    |
    |     N X |    |
    |   s|ubstate  |
    |    |  X N    |
    |    |   /     |
    ----N|- M- N----
         |/
         N

    :param cage: A Cage object
    :param substrate: A Substrate object
    :return: List of xyzs with both cage and substrate xyzs
    """
    logger.info('Adding substrate with 2(?) heteroatoms aligned along the M–M vector')

    cage_coords = xyz2coord(cage.xyzs)
    subst_coords = xyz2coord(substrate.xyzs)
    subst_x_id_0, subst_x_id_1 = substrate.x_x_atom_ids[0], substrate.x_x_atom_ids[1]
    substrate_center = calc_midpoint(subst_coords[subst_x_id_1], subst_coords[subst_x_id_0])

    # Translate the substrate to the M-M midpoints

    m_m_midpoint = calc_midpoint(cage_coords[cage.m_ids[1]], cage_coords[cage.m_ids[0]])
    cage_coords = [coord - m_m_midpoint for coord in cage_coords]  # Centre the cage on the M-M midpoint
    subst_coords = [coord - (m_m_midpoint + substrate_center) for coord in subst_coords]

    # Rotate the substrate so the X–X vector of the substrate is aligned with the M-M vector

    substrate_centre_x_vector = calc_normalised_vector(m_m_midpoint, subst_coords[subst_x_id_0])
    cage_m_center_vector = calc_normalised_vector(m_m_midpoint, cage_coords[cage.m_ids[1]])

    theta = np.arccos(np.dot(substrate_centre_x_vector, cage_m_center_vector))
    normal = np.cross(substrate_centre_x_vector, cage_m_center_vector)
    rot_matrix = rotation_matrix(normal, theta)
    subst_coords = [np.matmul(rot_matrix, coord) for coord in subst_coords]

    rot_substrate_coords = rot_minimise_repulsion(subst_coords, cage_coords, rot_axes=[cage_m_center_vector])
    xyzs = cat_cage_subst_coords(cage, substrate, cage_coords, rot_substrate_coords)

    return xyzs


def add_substrate_x(cage, substrate, m_x_dist=3.5):
    """
    Add a substrate with at least one heteroatom (X) which can form C(2)H hydrogen bonds, pointing at a M atom
    will just pick the first reasonable conformer...

    :param cage: A cage object
    :param substrate: A substrate object
    :param m_x_dist: Guess at the M-X bond distance in Å
    :return:
    """
    logger.info('Adding substrate with 1(?) heteroatom aligned along the M–M vector')

    cage_coords = xyz2coord(cage.xyzs)
    subst_coords = xyz2coord(substrate.xyzs)

    for x_atom_id in substrate.x_atom_ids:

        # Shift the M and X atoms to (0.0, 0.0, 0.0), shift the X atom to the first M atom
        cage_coords = [coord - cage_coords[cage.m_ids[0]] for coord in cage_coords]
        subst_coords = [coord - subst_coords[x_atom_id] for coord in subst_coords]

        # Assume the X atom is monovalent.... and bonded to an element E
        e_atom_id = get_closest_bonded_atom_id(substrate.xyzs, x_atom_id)
        e_x_bond_vector = calc_normalised_vector(subst_coords[e_atom_id], subst_coords[x_atom_id])
        m_m_vector = calc_normalised_vector(cage_coords[cage.m_ids[1]], cage_coords[cage.m_ids[0]])

        # Rotate the substrate to the E-X bond is inline with the M-M vector
        theta = np.arccos(np.dot(e_x_bond_vector, m_m_vector))
        normal = np.cross(e_x_bond_vector, m_m_vector)
        rot_matrix = rotation_matrix(normal, theta)
        subst_coords = [np.matmul(rot_matrix, coord) for coord in subst_coords]

        # Translate the substrate so the M-X distance is as defined
        subst_coords = [coord - (m_x_dist * m_m_vector) for coord in subst_coords]

        # Rotate in the M-M axis to minimise the steric repulsion
        subst_coords = rot_minimise_repulsion(subst_coords, cage_coords, rot_axes=[m_m_vector])

        xyzs = cat_cage_subst_coords(cage, substrate, cage_coords, subst_coords)

        # TODO loop over possible binding atoms, perhaps printing all possible binding modes

        return xyzs


def add_substrate_com(cage, substrate):
    """
    Add a substrate the centre of a cage defined by its centre of mass (com)
    :param cage: A Cage object
    :param substrate: A Substrate object
    :return:
    """
    logger.info('Adding substrate with 0(?) heteroatoms to the cage COM')

    cage_coords = xyz2coord(cage.xyzs)
    subst_coords = xyz2coord(substrate.xyzs)

    subst_com = calc_com(substrate.xyzs)

    m_m_midpoint = calc_midpoint(cage_coords[cage.m_ids[1]], cage_coords[cage.m_ids[0]])
    subst_coords = [coord - (m_m_midpoint + subst_com) for coord in subst_coords]

    rot1 = calc_normalised_vector(m_m_midpoint + subst_com, cage_coords[cage.m_ids[0]])
    rot2 = calc_normalised_vector(rot1[0] * rot1, np.array([1.0, 0.0, 0.0]))
    rot3 = calc_normalised_vector(rot1[1] * rot1 + rot2[1] * rot2, np.array([0.0, 1.0, 0.0]))

    subst_coords = rot_minimise_repulsion(subst_coords, cage_coords, rot_axes=[rot1, rot2, rot3], n_rot_steps=100)

    xyzs = cat_cage_subst_coords(cage, substrate, cage_coords, subst_coords)

    return xyzs


class Template:

    def set_geom(self, r):
        m_m_vec = self.Metals.coords[0] - self.Metals.coords[1]
        m_m_dist = np.linalg.norm(m_m_vec)
        m_m_normal_vec = m_m_vec / m_m_dist

        # Shift the second metal to a distance 2*r from the first
        self.Metals.coords[1] += (2.0 * r - m_m_dist) * m_m_normal_vec
        for i, ligand in enumerate(self.Ligands):
            for j in range(3, 6):
                self.Ligands[i].x_coords[j] += (2.0 * r - m_m_dist) * m_m_normal_vec
                self.Ligands[i].x_xyzs[j] = [self.Ligands[i].x_xyzs[j][0]] + self.Ligands[i].x_coords[j].tolist()

        self.Metals.xyzs = [[self.Metals.xyzs[i][0]] + self.Metals.coords[i].tolist() for i in range(2)]

    def __init__(self, r=None):
        self.Metals = Metals()
        self.Ligands = [Ligand1(), Ligand2(), Ligand3(), Ligand4()]

        if r is not None:
            self.set_geom(r)


class Metals:
    def __init__(self):
        self.xyzs = [['M', -6.05381, 0.00096, -0.00072],
                     ['M', 6.05381, -0.00054, 0.00066]]
        self.coords = xyz2coord(self.xyzs)


class Ligand1:

    def __init__(self):
        exe1_xyzs = [['E',  -7.18579, 1.92330, -1.92897],
                     ['X',  -6.03691, 1.43490, -1.43693],
                     ['E',  -4.86407, 1.88909, -1.88902]]

        exe2_xyzs = [['E',  7.18575, 1.92345, -1.92582],
                     ['X',  6.03680, 1.43433, -1.43465],
                     ['E',  4.86402, 1.88844, -1.88696]]

        self.x_xyzs = exe1_xyzs + exe2_xyzs
        self.x_coords = xyz2coord(self.x_xyzs)


class Ligand2:

    def __init__(self):
        exe1_xyzs = [['E', -7.18452, -1.92409, -1.92674],
                     ['X', -6.03583, -1.43469, -1.43525],
                     ['E', -4.86281, -1.88906, -1.88667]]

        exe2_xyzs = [['E', 7.18706, -1.92346, -1.92640],
                     ['X', 6.03792, -1.43522, -1.43481],
                     ['E', 4.86532, -1.89002, -1.88691]]

        self.x_xyzs = exe1_xyzs + exe2_xyzs
        self.x_coords = xyz2coord(self.x_xyzs)


class Ligand3:

    def __init__(self):
        exe1_xyzs = [['E', -7.18570, -1.92117, 1.92765],
                     ['X', -6.03676, -1.43298, 1.43551],
                     ['E', -4.86397, -1.88744, 1.88745]]

        exe2_xyzs = [['E', 7.18588, -1.92427, 1.92748],
                     ['X', 6.03698, -1.43539, 1.43596],
                     ['E', 4.86416, -1.88966, 1.88800]]

        self.x_xyzs = exe1_xyzs + exe2_xyzs
        self.x_coords = xyz2coord(self.x_xyzs)


class Ligand4:

    def __init__(self):
        exe1_xyzs = [['E', -7.18713, 1.92581, 1.92437],
                     ['X', -6.03797, 1.43657, 1.43382],
                     ['E', -4.86539, 1.89096, 1.88637]]

        exe2_xyzs = [['E',  7.18444, 1.92318, 1.92804],
                     ['X',  6.03577, 1.43419, 1.43609],
                     ['E',  4.86273, 1.88841, 1.88764]]

        self.x_xyzs = exe1_xyzs + exe2_xyzs
        self.x_coords = xyz2coord(self.x_xyzs)
