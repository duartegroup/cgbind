import numpy as np
from cgbind.log import logger
from cgbind.config import Config
from cgbind.atoms import heteroatoms
from cgbind.geom import calc_distance_matrix
from cgbind.geom import get_neighbour
from cgbind.geom import xyz2coord
from cgbind.geom import calc_normalised_vector as calc_norm_vec
from cgbind.build import get_cost_xyzs
from scipy.optimize import minimize


def build(cage, linker):
    """
    Build a M4L6 metallocage


             M
           `-.`
          `-.``.`
         .--``````
        .--.``````.`
      `----```````````
     `-----````````````
    .-----.```o..`r````
   M------````-..----...``
    .----.``````````...----M
     .---``````````.```````
      .-.``.``````
       ``M

    :param cage: A cage object
    :param linker: A linker object
    """
    logger.info('Building a M4L6 metallocage {}'.format(cage.name))

    opt_res = minimize(get_cost_xyzs, x0=np.array([1.0]), method='BFGS',
                       args=(cage, linker, Template(), get_xeex_coords(linker), False))

    cost, xyzs = get_cost_xyzs(r=opt_res.x[0], cage=cage, linker=linker, template=Template(),
                               coords_to_fit=get_xeex_coords(linker), return_xyzs=True)

    logger.info('Final centre-M dist {} with cost value {}'.format(np.round(opt_res.x[0], 1),
                                                                   np.round(cost, 1)))

    if len(xyzs) != cage.arch.n_metals + cage.arch.n_linkers * linker.n_atoms:
        logger.error('Could not generate M4L6 metallocage')
    else:
        cage.xyzs = xyzs

    return None


def has_correct_xee_angle(xyzs, x, e1, e2, correct_angle_deg=120, tolerance=0.2):
    """
    For a reasonable xeex motif the xee angles should be ~ 120 degrees

    X1        X2
     \       /
      E1----E2
    :param xyzs:
    :param x:
    :param e1:
    :param e2:
    :param correct_angle_deg:
    :param tolerance: Relative tolerance on the correct angle
    :return:
    """

    e1_x_vec = calc_norm_vec(xyz2coord(xyzs[e1]), xyz2coord(xyzs[x]))
    e1_e2_vec = calc_norm_vec(xyz2coord(xyzs[e1]), xyz2coord(xyzs[e2]))
    dot_vecs = np.dot(e1_x_vec, e1_e2_vec)

    if dot_vecs > 1.0:                              # Save numpy.arccos from floating point precision
        dot_vecs = 1.0
    angle_deg = np.rad2deg(np.arccos(dot_vecs))
    return correct_angle_deg * (1.0 - tolerance) < angle_deg < correct_angle_deg * (1.0 + tolerance)


def is_xeex_motif(xyzs, x1, e1, e2, x2, max_xx_dist, max_ee_dist):
    """
    Is a set of atom ids a xeex motif, which will be true if both x1 and x2 are heteroatoms, they are separated by a
    small enough distance and the e1 and e2 atoms are bonded (separated by a small enough distance)

    :param xyzs: List of xyzs
    :param x1: atom id
    :param e1: atom id
    :param e2: atom id
    :param x2: atom id
    :param max_xx_dist: Maximum distance for the heteroatoms to be from each other (Å)
    :param max_ee_dist: Maximum distance for the adjoining atoms to be from each other (Å)
    :return:
    """

    x_x_dist = np.linalg.norm(xyz2coord(xyzs[x1]) - xyz2coord(xyzs[x2]))

    if xyzs[x1][0] in heteroatoms and xyzs[x2][0] in heteroatoms:
        if x_x_dist < max_xx_dist:
            if np.linalg.norm(xyz2coord(xyzs[e1]) - xyz2coord(xyzs[e2])) < max_ee_dist:
                if e1 != e2:
                    if has_correct_xee_angle(xyzs, x1, e1, e2) and has_correct_xee_angle(xyzs, x2, e2, e1):
                        return True

    return False


def get_confomers_with_flat_xeex_motifs(conf_xyzs, conf_ids, xeex_motifs, tol=0.1):
    """

    :param conf_xyzs: (list(list)) list of xyzs
    :param conf_ids: (list) list of conformers ids
    :param xeex_motifs: (list(list))
    :param tol: (float) relative tolerance on the ideal angle between normal vectors = 0º
    :return:
    """
    ideal_abs_cos_angle_norm_vec = 1.0

    flat_conf_ids = []
    x1_1, e1_1, e2_1, x2_1 = xeex_motifs[0]
    x1_2, e1_2, e2_2, x2_2 = xeex_motifs[1]

    for conf_id in conf_ids:
        xyzs = conf_xyzs[conf_id]
        e11_x11_vec = calc_norm_vec(xyz2coord(xyzs[e1_1]), xyz2coord(xyzs[x1_1]))
        e11_e21_vec = calc_norm_vec(xyz2coord(xyzs[e1_1]), xyz2coord(xyzs[e2_1]))
        norm1 = np.cross(e11_x11_vec, e11_e21_vec) / np.linalg.norm(np.cross(e11_x11_vec, e11_e21_vec))

        e12_x12_vec = calc_norm_vec(xyz2coord(xyzs[e1_2]), xyz2coord(xyzs[x1_2]))
        e12_e22_vec = calc_norm_vec(xyz2coord(xyzs[e1_2]), xyz2coord(xyzs[e2_2]))
        norm2 = np.cross(e12_x12_vec, e12_e22_vec) / np.linalg.norm(np.cross(e12_x12_vec, e12_e22_vec))

        cos_angle = np.dot(norm1, norm2)
        if (1.0 - tol) * ideal_abs_cos_angle_norm_vec < np.abs(cos_angle) < (1.0 + tol) * ideal_abs_cos_angle_norm_vec:
            flat_conf_ids.append(conf_id)

    return flat_conf_ids


def get_and_set_xeex_motifs(linker, xyzs, max_x_x_dist=4.0, max_e_e_dist=2.0):

    linker.set_com()
    n_xeex_motifs = 0
    xeex_atom_ids = []

    distance_matrix = calc_distance_matrix(xyzs)
    for i in range(linker.n_atoms):
        for j in range(linker.n_atoms):
            if i > j:
                xi_nn, xi_nnn = get_neighbour(i, distance_matrix, 1), get_neighbour(i, distance_matrix, 2)
                xj_nn, xj_nnn = get_neighbour(j, distance_matrix, 1), get_neighbour(j, distance_matrix, 2)

                for e1 in [xi_nn, xi_nnn]:
                    for e2 in [xj_nn, xj_nnn]:
                        if is_xeex_motif(xyzs, i, e1, e2, j, max_x_x_dist, max_e_e_dist):
                            n_xeex_motifs += 1
                            # The XEEX motifs are defined from the linker COM to the furthest atom
                            if (np.linalg.norm(xyz2coord(xyzs[j]) - linker.com) >
                                    np.linalg.norm(xyz2coord(xyzs[i]) - linker.com)):
                                xeex_atom_ids.append([i, e1, e2, j])
                            else:
                                xeex_atom_ids.append([j, e2, e1, i])

    if n_xeex_motifs != 2:
        logger.critical('Cannot handle >2 XEEX motifs in the linker yet')

    linker.xeex_motifs = xeex_atom_ids
    return xeex_atom_ids


def get_best_linker_conformer(linker):
    """

    For a linker to be suitable for constructing a tetrahedral metallocage of the form M4L6

     M ------ M
     |\      /|
     |  \  /  |                  (<- that is meant to be a tetrahedron)
     |   /\   |
     | /    \ |
     M ------ M

     The metal atoms are assumed to be octahedrally coordinated by bidentate linker termini, with structure


      X     X
      \____/
        |___________
                  __|__
                 /     \
                X      X


    This function with first search for X-E-E-X motifs where the X-E-E angles are 120 degrees, so that the metal
    has 90 degree angles when octahedrally coordinated,

    The best confomer will have vectors (vx) which are roughly inverses with cos(theta) = -1

         vx
         ^
      X  |   X
      \__|__/


    :param linker: A linker object
    :param max_x_x_dist:
    :param max_e_e_dist:
    :return: Best conformer id
    """
    best_conf_id = None

    xeex_atom_ids = get_and_set_xeex_motifs(linker, xyzs=linker.conf_xyzs[0])
    flat_conf_ids = get_confomers_with_flat_xeex_motifs(linker.conf_xyzs, linker.conf_ids, xeex_atom_ids)

    curr_cos_theta = 0.0
    ideal_cos_theta = -1.0                      #
    for conf_id in flat_conf_ids:
        xyzs = linker.conf_xyzs[conf_id]

        x1_e1_vec1 = xyz2coord(xyzs[xeex_atom_ids[0][0]]) - xyz2coord(xyzs[xeex_atom_ids[0][1]])
        x2_e2_vec1 = xyz2coord(xyzs[xeex_atom_ids[0][3]]) - xyz2coord(xyzs[xeex_atom_ids[0][2]])
        xe_xe_ang1 = np.arccos(np.dot(x1_e1_vec1, x2_e2_vec1) / (np.linalg.norm(x1_e1_vec1)*np.linalg.norm(x2_e2_vec1)))
        v1 = (x1_e1_vec1 + x2_e2_vec1) / np.linalg.norm(x1_e1_vec1 + x2_e2_vec1)

        x1_e1_vec2 = xyz2coord(xyzs[xeex_atom_ids[1][0]]) - xyz2coord(xyzs[xeex_atom_ids[1][1]])
        x2_e2_vec2 = xyz2coord(xyzs[xeex_atom_ids[1][3]]) - xyz2coord(xyzs[xeex_atom_ids[1][2]])
        xe_xe_ang2 = np.arccos(np.dot(x1_e1_vec2, x2_e2_vec2) / (np.linalg.norm(x1_e1_vec2)*np.linalg.norm(x2_e2_vec2)))
        v2 = (x1_e1_vec2 + x2_e2_vec2) / np.linalg.norm(x1_e1_vec2 + x2_e2_vec2)

        cos_theta = np.dot(v1, v2)

        if xe_xe_ang1 < np.pi/2.0 and xe_xe_ang2 < np.pi/2.0:
            if np.abs(cos_theta - ideal_cos_theta) < np.abs(curr_cos_theta - ideal_cos_theta):
                best_conf_id = conf_id
                curr_cos_theta = cos_theta

    if not Config.suppress_print:
        print("{:<30s}{:<50s}{:>10s}".format('Best conformer of', linker.name, 'Found'))

    return best_conf_id


def get_xeex_coords(linker):
    """
    Get the coordinates of the twp XEEX motifs in the tmp_linker in order as they appear from one end of the tmp_linker
    to the other
    :param linker: (object)
    :return:
    """
    return np.array([xyz2coord(linker.xyzs[i]) for motif in linker.xeex_motifs for i in motif])


class Template:

    def set_geom(self, r, template_atoms_cuttoff=3.0):
        centroid = np.average(xyz2coord(self.Metals.xyzs), axis=0)

        for i, metal_coord in enumerate(self.Metals.coords):
            # Define the current centroid – metal vector a
            c_m_vec = calc_norm_vec(coord1=metal_coord, coord2=centroid)
            curr_c_m_dist = np.linalg.norm(centroid - metal_coord)

            # Distance to translate along c_m_vec
            trns = (r - curr_c_m_dist)

            for j, ligand in enumerate(self.Ligands):
                for k, l_coord in enumerate(ligand.x_coords):
                    l_atom_m_dist = np.linalg.norm(metal_coord - l_coord)
                    if l_atom_m_dist < template_atoms_cuttoff:
                        self.Ligands[j].x_coords[k] += trns * c_m_vec
                        self.Ligands[j].x_xyzs[k] = ([self.Ligands[j].x_xyzs[k][0]] +
                                                     self.Ligands[j].x_coords[k].tolist())

            self.Metals.coords[i] += trns * c_m_vec
            self.Metals.xyzs[i] = [self.Metals.xyzs[i][0]] + self.Metals.coords[i].tolist()

    def __init__(self, r=None):
        self.Metals = Metals()
        self.Ligands = [Ligand1(), Ligand2(), Ligand3(), Ligand4(), Ligand5(), Ligand6()]

        if r is not None:
            self.set_geom(r)


class Metals:
    def __init__(self):
        self.xyzs = [['M', 0.00000, 0.00000, 12.69480],
                     ['M', 3.83600, 6.33800, 23.13010],
                     ['M', 3.57090, -6.49110, 23.13010],
                     ['M', -7.40690, 0.15310, 23.13010]]
        self.coords = xyz2coord(self.xyzs)


class Ligand1:

    def __init__(self):
        xeex1_xyzs = [['X', -0.11360, 1.59760, 11.59410],  # 1,1: X1
                      ['E', -0.97970, 2.52000, 12.04870],  # 1,1: E1
                      ['E', -1.62480, 2.25460, 13.26720],  # 1,1: E2
                      ['X', -1.21420, 1.11020, 13.84450]]  # 1,1: X2

        xeex2_xyzs = [['X', -6.72220, -1.15890, 21.78490],  # 1,2: X1
                      ['E', -7.26930, -2.38460, 21.88140],  # 1,2: E1
                      ['E', -8.14270, -2.49650, 22.96240],  # 1,2: E2
                      ['X', -8.41060, -1.39720, 23.69840]]  # 1,2: X2

        self.x_xyzs = xeex1_xyzs + xeex2_xyzs
        self.x_coords = xyz2coord(self.x_xyzs)


class Ligand2:
    def __init__(self):
        xeex1_xyzs = [['X', -1.32670, -0.89720, 11.59410],  # 2,1: X1
                      ['E', -1.69250, -2.10840, 12.04870],  # 2,1: E1
                      ['E', -1.14020, -2.53440, 13.26720],  # 2,1: E2
                      ['X', -0.35430, -1.60660, 13.84450]]  # 2,1: X2

        xeex2_xyzs = [['X', 5.41530, -6.58520, 23.69840],  # 2,2: X1
                      ['E', 6.23340, -5.80350, 22.96240],  # 2,2: E1
                      ['E', 5.69980, -5.10310, 21.88140],  # 2,2: E2
                      ['X', 4.36470, -5.24210, 21.78490]]  # 2,2: X2
        self.x_xyzs = xeex1_xyzs + xeex2_xyzs
        self.x_coords = xyz2coord(self.x_xyzs)


class Ligand3:
    def __init__(self):
        xeex1_xyzs = [['X', 1.44030, -0.70040, 11.59410],  # 3,1: X1
                      ['E', 2.67220, -0.41160, 12.04870],  # 3,1: E1
                      ['E', 2.76500, 0.27980, 13.26720],  # 3,1: E2
                      ['X', 1.56850, 0.49640, 13.84450]]  # 3,1: X2

        xeex2_xyzs = [['X', 2.99530, 7.98230, 23.69840],  # 3,2: X1
                      ['E', 1.90930, 8.30000, 22.96240],  # 3,2: E1
                      ['E', 1.56960, 7.48770, 21.88140],  # 3,2: E2
                      ['X', 2.35750, 6.40100, 21.78490]]  # 3,2: X2
        self.x_xyzs = xeex1_xyzs + xeex2_xyzs
        self.x_coords = xyz2coord(self.x_xyzs)


class Ligand4:
    def __init__(self):
        xeex1_xyzs = [['X', -7.70871, 1.16494, 24.58487],  # 4,1: X1
                      ['E', -6.96610, 1.01090, 25.71820],  # 4,1: E1
                      ['E', -5.83940, 0.19320, 25.47620],  # 4,1: E2
                      ['X', -5.80510, -0.21660, 24.1760]]  # 4,1: X2

        xeex2_xyzs = [['X', 3.74360, -8.11950, 22.17570],  # 4,2: X1
                      ['E', 2.63570, -8.50760, 21.48390],  # 4,2: E1
                      ['E', 1.51740, -7.66100, 21.56240],  # 4,2: E2
                      ['X', 1.80090, -6.53460, 22.25590]]  # 4,2: X2
        self.x_xyzs = xeex1_xyzs + xeex2_xyzs
        self.x_coords = xyz2coord(self.x_xyzs)


class Ligand5:
    def __init__(self):
        xeex1_xyzs = [['X', 2.99220, -7.34870, 24.79090],  # 5,1: X1
                      ['E', 2.60760, -6.53820, 25.71820],  # 5,1: E1
                      ['E', 2.75250, -5.15370, 25.47620],  # 5,1: E2
                      ['X', 3.09010, -4.91900, 24.17600]]  # 5,1: X2

        xeex2_xyzs = [['X', 5.15990, 7.30180, 22.17570],  # 5,2: X1
                      ['E', 6.05000, 6.53640, 21.48390],  # 5,2: E1
                      ['E', 5.87590, 5.14470, 21.56240],  # 5,2: E2
                      ['X', 4.75870, 4.82700, 22.25590]]  # 5,2: X2
        self.x_xyzs = xeex1_xyzs + xeex2_xyzs
        self.x_coords = xyz2coord(self.x_xyzs)


class Ligand6:
    def __init__(self):
        xeex1_xyzs = [['X', 4.86810, 6.26570, 24.79090],  # 6,1: X1
                      ['E', 4.35850, 5.52730, 25.71820],  # 6,1: E1
                      ['E', 3.08700, 4.96050, 25.47620],  # 6,1: E2
                      ['X', 2.71490, 5.13560, 24.17600]]  # 6,1: X2

        xeex2_xyzs = [['X', -8.90350, 0.81770, 22.17570],  # 6,2: X1
                      ['E', -8.68570, 1.97120, 21.48390],  # 6,2: E1
                      ['E', -7.39340, 2.51640, 21.56240],  # 6,2: E2
                      ['X', -6.55960, 1.70770, 22.25590]]  # 6,2: X2
        self.x_xyzs = xeex1_xyzs + xeex2_xyzs
        self.x_coords = xyz2coord(self.x_xyzs)
