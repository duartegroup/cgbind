import sys
import numpy as np
from rdkit import Chem
from rdkit.Numerics import rdAlignment
from rdkit.Chem import AllChem
from copy import deepcopy
from .config import Config
from .atoms import heteroatoms
from .geom import calc_distance_matrix
from .geom import get_neighbour
from .geom import xyz2coord
from .geom import calc_normalised_vector as calc_norm_vec
from .geom import calc_midpoint
from .geom import rotation_matrix
from .geom import molblock2xyzs
from .geom import is_geom_reasonable


def build(cage, linker):
    """
    Build a M4L6 metallocage

    :param cage: A cage object
    :param linker: A linker object
    :return:
    """

    n_metals, n_ligands = 4, 6
    xyzs = []

    template_center = np.array([np.average([xyz2coord(xyz)[i] for xyz in Template.Metals.xyzs]) for i in range(3)])
    template_radius = np.average([np.linalg.norm(xyz2coord(xyz) - template_center) for xyz in Template.Metals.xyzs])
    center_metal_vectors = [xyz2coord(xyz) - template_center for xyz in Template.Metals.xyzs]

    linker_xeex_xyzs = [linker.xyzs[linker.xeex_motifs[i][j]] for i in range(2) for j in range(4)]
    template_linker_xyzs = Template.Ligands[0].xeex1_xyzs + Template.Ligands[0].xeex2_xyzs

    rmsd, transform_mat = rdAlignment.GetAlignmentTransform(xyz2coord(template_linker_xyzs).tolist(),
                                            xyz2coord(linker_xeex_xyzs).tolist())
    tmp_linker_obj = deepcopy(linker.mol_obj)
    AllChem.TransformMol(tmp_linker_obj, transform_mat)

    transformed_coords = [[tmp_linker_obj.GetConformer().GetAtomPosition(i)[j] for j in range(3)] for i in range(linker.n_atoms)]
    transformed_xyzs = [[linker.xyzs[i][0]] + transformed_coords[i] for i in range(linker.n_atoms)]

    from .input_output import xyzs2xyzfile
    xyzs2xyzfile(Template.Metals.xyzs + transformed_xyzs, basename='tmp2')

    exit()

    return xyzs


def uff_cage_opt(xyzs, cage, linker, n_ligands=6, n_metals=4, tmp_m_atom='S'):
    """
    Do a universal forcefield (UFF) optimisation of a temporary M4L6 cage by converting it first to an RDKit mol object,
    inhering the temporary 'metal' atoms, adding bonds, modifying the formal charges and then attempting the
    optimisation

    :param xyzs:
    :param cage:
    :param linker:
    :param n_ligands:
    :param n_metals:
    :param tmp_m_atom:
    :return:
    """

    tmp_m_atom_obj = Chem.MolFromSmiles(tmp_m_atom)
    tmp_cage_rdkit_obj = combine_rdkit_objs(linker.mol_obj, tmp_m_atom_obj)
    cage_rdkit_molblock_lines = Chem.MolToMolBlock(tmp_cage_rdkit_obj).split('\n')
    new_molblock = insert_new_coords_mol_block(cage_rdkit_molblock_lines, xyzs)
    tmp_cage_rdkit_obj = Chem.MolFromMolBlock(new_molblock, removeHs=False)
    tmp_cage_dist_matrix = calc_distance_matrix(molblock2xyzs(Chem.MolToMolBlock(tmp_cage_rdkit_obj)))
    bonds_to_add = gen_bonds_to_add(xyzs, linker, n_metals, n_ligands, tmp_cage_dist_matrix)
    tmp_cage_rdkit_obje = Chem.EditableMol(tmp_cage_rdkit_obj)
    [tmp_cage_rdkit_obje.AddBond(*bond, order=Chem.rdchem.BondType.SINGLE) for bond in bonds_to_add]
    tmp_cage_rdkit_obj = tmp_cage_rdkit_obje.GetMol()
    add_formal_charges(tmp_cage_rdkit_obj)
    tmp_cage_rdkit_obj.UpdatePropertyCache()
    Chem.GetSymmSSSR(tmp_cage_rdkit_obj)
    tmp_cage_rdkit_obj.GetRingInfo().NumRings()

    try:
        AllChem.UFFOptimizeMolecule(tmp_cage_rdkit_obj)
        tmp_cage_rdkit_obj = Chem.AddHs(tmp_cage_rdkit_obj)
    except RuntimeError or ValueError:
        return None

    xyzs = molblock2xyzs(Chem.MolToMolBlock(tmp_cage_rdkit_obj))

    # Replace the tmp_m_atom with the actual metal atoms in the list of xys
    for i in range(n_metals):
        xyzs[-(1+i)][0] = cage.metal

    return xyzs


def has_ok_m_x_distances(xyzs, linker, n_ligands=6, n_metals=4, max_mx_dist=4.0):

    dist_matrix = calc_distance_matrix(xyzs)
    bonds_to_add = gen_bonds_to_add(xyzs, linker, n_metals, n_ligands, dist_matrix)
    bond_dists_to_add = [dist_matrix[bond[0], bond[1]] for bond in bonds_to_add]

    if all([dist < max_mx_dist for dist in bond_dists_to_add]):
        return True
    else:
        return False


def add_linker_xyzs(xyzs, linker, first_linker_coords, axes, thetax=0.0, theray=0.0, thetaz=0.0, thetaz2=None,
                    trnsx=0.0, trnsy=0.0, trnsz=0.0):

    x_ax, y_ax, z_ax = axes

    linkerx_coords = deepcopy(first_linker_coords)
    linkerx_coords = [np.matmul(rotation_matrix(z_ax, thetaz), coord) for coord in linkerx_coords]
    linkerx_coords = [np.matmul(rotation_matrix(y_ax, theray), coord) for coord in linkerx_coords]
    linkerx_coords = [np.matmul(rotation_matrix(x_ax, thetax), coord) for coord in linkerx_coords]
    if thetaz2:
        linkerx_coords = [np.matmul(rotation_matrix(z_ax, thetaz2), coord) for coord in linkerx_coords]
    linkerx_coords = [coord + z_ax * trnsz for coord in linkerx_coords]
    linkerx_coords = [coord + y_ax * trnsy for coord in linkerx_coords]
    linkerx_coords = [coord + x_ax * trnsx for coord in linkerx_coords]

    xyzs += [[linker.xyzs[i][0]] + linkerx_coords[i].tolist() for i in range(linker.n_atoms)]

    return xyzs


def add_metal_xyzs(xyzs, cage, linker, l1_coords, l2_coords):

    for l_coords in [l1_coords, l2_coords]:
        for i in range(2):

            xeex_motif = linker.xeex_motifs[i]
            e11_x11_vec = calc_norm_vec(l_coords[xeex_motif[1]], l_coords[xeex_motif[0]])
            e21_x21_vec = calc_norm_vec(l_coords[xeex_motif[2]], l_coords[xeex_motif[3]])
            sum_xe_vec = e11_x11_vec + e21_x21_vec
            metal_coord = calc_midpoint(l_coords[xeex_motif[0]], l_coords[xeex_motif[3]]) + sum_xe_vec

            xyzs += [[cage.metal] + metal_coord.tolist()]

    return xyzs


def combine_rdkit_objs(linker_obj, metal_obj):
    """
    rdkit can only combine objects in pairs, so loop through all the metal and linker objects and combine
    :param linker_obj:
    :param metal_obj:
    :return:
    """
    n_metals, n_ligands = 4, 6

    rdkit_obj = Chem.CombineMols(linker_obj, linker_obj)
    for i in range(n_ligands - 2):
        rdkit_obj = Chem.CombineMols(rdkit_obj, linker_obj)

    for i in range(n_metals):
        rdkit_obj = Chem.CombineMols(rdkit_obj, metal_obj)

    return rdkit_obj


def insert_new_coords_mol_block(molblock_lines, xyzs):

    new_molblock_lines = []
    n_line = 0
    for line in molblock_lines:
        if len(line.split()) == 16 and line.split()[0][-1].isdigit():
            x, y, z = xyzs[n_line][1:]
            x, y, z = str(np.round(x, 4)), str(np.round(y, 4)), str(np.round(z, 4))
            new_line = ('{:>10s}{:>10s}{:>10s}{:>3s}{:>3s}{:>3s}{:>3s}{:>3s}{:>3s}{:>3s}'
                        '{:>3s}{:>3s}{:>3s}{:>3s}{:>3s}{:>3s}').format(x, y, z, *line.split()[3:])

            new_molblock_lines.append(new_line)
            n_line += 1
        else:
            new_molblock_lines.append(line)
    new_molblock = '\n'.join(new_molblock_lines)

    return new_molblock


def gen_axes(linker, l_coords):
    """
    Generate the x, y, z axes centred on the middle of a linker, with the z axis normal to the xeex motifs, then the
    x a y axes orthogonal to those
    :param linker: Linker object
    :param l_coords: Linker coordinates (numpy array of [x, y, z] values for each atom)
    :return:
    """

    e11_x11_vec = calc_norm_vec(l_coords[linker.xeex_motifs[0][1]], l_coords[linker.xeex_motifs[0][0]])
    e11_e21_vec = calc_norm_vec(l_coords[linker.xeex_motifs[0][1]], l_coords[linker.xeex_motifs[0][2]])
    e21_x21_vec = calc_norm_vec(l_coords[linker.xeex_motifs[0][2]], l_coords[linker.xeex_motifs[0][3]])

    z_axis = np.cross(e11_x11_vec, e11_e21_vec) / np.linalg.norm(np.cross(e11_x11_vec, e11_e21_vec))
    sum_xe_vec = e11_x11_vec + e21_x21_vec
    x_axis = calc_norm_vec(np.dot(sum_xe_vec, z_axis) * z_axis, sum_xe_vec)                 # Gram–Schmidt
    guess_yax = np.array([1.0, 0.0, 0.0])
    y_axis = calc_norm_vec(np.dot(guess_yax, z_axis) * z_axis + np.dot(guess_yax, x_axis) * x_axis, guess_yax)

    return [x_axis, y_axis, z_axis]


def gen_bonds_to_add(xyzs, linker, n_metals, n_ligands, dist_matrix):
    """
    For a temporary cage structure add the bonds between the 'M' atoms and the surrounding xeex motifs, so the
    molecular mechanics optimiser has an idea of what to do...
    :param xyzs:
    :param linker:
    :param n_metals:
    :param n_ligands:
    :param dist_matrix:
    :return:
    """

    bonds_to_add = []
    for tmp_m_id in range(len(xyzs) - n_metals, len(xyzs)):

        closest_x_ids = []
        for i in range(1, len(xyzs)):
            if len(closest_x_ids) == n_ligands:
                break
            neighbour_id = get_neighbour(tmp_m_id, dist_matrix, i)
            for xeex_motif in linker.xeex_motifs:
                if neighbour_id in [xeex_motif[j] + (n * linker.n_atoms) for j in [0, 3] for n in range(n_ligands)]:

                    closest_x_ids.append(neighbour_id)
                    bonds_to_add.append([tmp_m_id, neighbour_id])

    return bonds_to_add


def add_formal_charges(mol_obj):
    """
    Alter the formal charges on atoms so MM optimiser can cope...
    :param mol_obj:
    :return:
    """

    mol_obj.UpdatePropertyCache(strict=False)
    for atom in mol_obj.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetExplicitValence() == 4 and atom.GetFormalCharge() == 0:
            atom.SetFormalCharge(1)
        if atom.GetAtomicNum() == 7 and atom.GetExplicitValence() == 3 and atom.GetFormalCharge() == -1:
            atom.SetFormalCharge(0)
        if atom.GetAtomicNum() == 8 and atom.GetExplicitValence() == 3 and atom.GetFormalCharge() == 0:
            atom.SetFormalCharge(1)
        if atom.GetAtomicNum() == 8 and atom.GetExplicitValence() == 2 and atom.GetFormalCharge() == -1:
            atom.SetFormalCharge(0)
    return 0


def calc_cube_side_length(coords, xeex_motifs):
    """
    Calculate the rough side length of the M4L6 cube
    :param coords: Linker coordinates (np array)
    :param xeex_motifs:
    :return:
    """

    p1 = calc_midpoint(coords[xeex_motifs[0][0]], coords[xeex_motifs[0][3]])
    p2 = calc_midpoint(coords[xeex_motifs[1][0]], coords[xeex_motifs[1][3]])

    return np.linalg.norm(p1 - p2)


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

    :param conf_xyzs:
    :param conf_ids:
    :param xeex_motifs:
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


def get_best_linker_conformer(linker, max_x_x_dist=4.0, max_e_e_dist=2.0):
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

    n_xeex_motifs = 0
    xeex_atom_ids = []

    first_xyzs = linker.conf_xyzs[0]
    distance_matrix = calc_distance_matrix(first_xyzs)
    for i in range(linker.n_atoms):
        for j in range(linker.n_atoms):
            if i > j:
                xi_nn, xi_nnn = get_neighbour(i, distance_matrix, 1), get_neighbour(i, distance_matrix, 2)
                xj_nn, xj_nnn = get_neighbour(j, distance_matrix, 1), get_neighbour(j, distance_matrix, 2)

                for e1 in [xi_nn, xi_nnn]:
                    for e2 in [xj_nn, xj_nnn]:
                        if is_xeex_motif(first_xyzs, i, e1, e2, j, max_x_x_dist, max_e_e_dist):
                            n_xeex_motifs += 1
                            xeex_atom_ids.append([i, e1, e2, j])

    # TODO allow for > 2 xeex motifs
    if n_xeex_motifs != 2:
        return None

    best_conf_id = 0
    linker.xeex_motifs = xeex_atom_ids

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


class Template(object):
    class Metals(object):
        xyzs = [['M', 0.00000, 0.00000, 12.69480],
                ['M', 3.83600, 6.33800, 23.13010],
                ['M', 3.57090, -6.49110, 23.13010],
                ['M', -7.40690, 0.15310, 23.13010]]

    class Ligand1(object):
        xeex1_xyzs = [['X', -0.11360, 1.59760, 11.59410],  # 1,1: X1
                      ['E', -0.97970, 2.52000, 12.04870],  # 1,1: E1
                      ['E', -1.62480, 2.25460, 13.26720],  # 1,1: E2
                      ['X', -1.21420, 1.11020, 13.84450]]  # 1,1: X2
        xeex2_xyzs = [['X', -6.72220, -1.15890, 21.78490],  # 1,2: X1
                      ['E', -7.26930, -2.38460, 21.88140],  # 1,2: E1
                      ['E', -8.14270, -2.49650, 22.96240],  # 1,2: E2
                      ['X', -8.41060, -1.39720, 23.69840]]  # 1,2: X2

    class Ligand2(object):
        xeex1_xyzs = [['X', -1.32670, -0.89720, 11.59410],  # 2,1: X1
                      ['E', -1.69250, -2.10840, 12.04870],  # 2,1: E1
                      ['E', -1.14020, -2.53440, 13.26720],  # 2,1: E2
                      ['X', -0.35430, -1.60660, 13.84450]]  # 2,1: X2
        xeex2_xyzs = [['X', 5.41530, -6.58520, 23.69840],  # 2,2: X1
                      ['X', 4.36470, -5.24210, 21.78490],  # 2,2: X2
                      ['E', 6.23340, -5.80350, 22.96240],  # 2,2: E1
                      ['E', 5.69980, -5.10310, 21.88140]]  # 2,2: E2

    class Ligand3(object):
        xeex1_xyzs = [['X', 1.44030, -0.70040, 11.59410],  # 3,1: X1
                      ['E', 2.67220, -0.41160, 12.04870],  # 3,1: E1
                      ['E', 2.76500, 0.27980, 13.26720],  # 3,1: E2
                      ['X', 1.56850, 0.49640, 13.84450]]  # 3,1: X2
        xeex2_xyzs = [['X', 2.99530, 7.98230, 23.69840],  # 3,2: X1
                      ['E', 1.90930, 8.30000, 22.96240],  # 3,2: E1
                      ['E', 1.56960, 7.48770, 21.88140],  # 3,2: E2
                      ['X', 2.35750, 6.40100, 21.78490]]  # 3,2: X2

    class Ligand4(object):
        xeex1_xyzs = [['X', -7.86030, 1.08310, 24.79090],  # 4,1: X1
                      ['E', -6.96610, 1.01090, 25.71820],  # 4,1: E1
                      ['E', -5.83940, 0.19320, 25.47620],  # 4,1: E2
                      ['X', -5.80510, -0.21660, 24.17600]]  # 4,1: X2

        xeex2_xyzs = [['X', 3.74360, -8.11950, 22.17570],  # 4,2: X1
                      ['E', 2.63570, -8.50760, 21.48390],  # 4,2: E1
                      ['E', 1.51740, -7.66100, 21.56240],  # 4,2: E2
                      ['X', 1.80090, -6.53460, 22.25590]]  # 4,2: X2

    class Ligand5(object):
        xeex1_xyzs = [['X', 2.99220, -7.34870, 24.79090],  # 5,1: X1
                      ['E', 2.60760, -6.53820, 25.71820],  # 5,1: E1
                      ['E', 2.75250, -5.15370, 25.47620],  # 5,1: E2
                      ['X', 3.09010, -4.91900, 24.17600]]  # 5,1: X2

        xeex2_xyzs = [['X', 5.15990, 7.30180, 22.17570],  # 5,2: X1
                      ['E', 6.05000, 6.53640, 21.48390],  # 5,2: E1
                      ['E', 5.87590, 5.14470, 21.56240],  # 5,2: E2
                      ['X', 4.75870, 4.82700, 22.25590]]  # 5,2: X2

    class Ligand6(object):
        xeex1_xyzs = [['X', 4.86810, 6.26570, 24.79090],  # 6,1: X1
                      ['E', 4.35850, 5.52730, 25.71820],  # 6,1: E1
                      ['E', 3.08700, 4.96050, 25.47620],  # 6,1: E2
                      ['X', 2.71490, 5.13560, 24.17600]]  # 6,1: X2

        xeex2_xyzs = [['X', -8.90350, 0.81770, 22.17570],  # 6,2: X1
                      ['E', -8.68570, 1.97120, 21.48390],  # 6,2: E1
                      ['E', -7.39340, 2.51640, 21.56240],  # 6,2: E2
                      ['X', -6.55960, 1.70770, 22.25590]]  # 6,2: X2

    Ligands = [Ligand1, Ligand2, Ligand3, Ligand4, Ligand5, Ligand6]


def old_build(cage, linker):
    """
    Build a M4L6 metallocage

    :param cage: A cage object
    :param linker: A linker object
    :return:
    """

    n_metals, n_ligands = 4, 6
    tmp_m_atom = 'S'            # Sulfur works nicely here, as it can support CN=6

    l1_coords = xyz2coord(linker.xyzs)
    cube_l = calc_cube_side_length(l1_coords, linker.xeex_motifs) * 1.2
    linker_midpoint = calc_midpoint(l1_coords[linker.xeex_motifs[0][0]], l1_coords[linker.xeex_motifs[1][0]])
    xax, yax, zax = gen_axes(linker, l1_coords)  # x, y, z axes

    # Shift the first linker to the origin (0.0, 0.0, 0.0)
    l1_coords = [coord - linker_midpoint for coord in l1_coords]

    # Add the xyz lines
    for axes in [[xax, yax, -zax], [-xax, yax, zax], [xax, -yax, zax], [-xax, -yax, zax], [xax, yax, -zax]]:

        l2_coords = [np.matmul(rotation_matrix(axes[0], np.pi), coord) for coord in deepcopy(l1_coords)]
        l2_coords = [coord + axes[2] * cube_l for coord in l2_coords]

        xyzs = []
        xyzs = add_linker_xyzs(xyzs, linker, l1_coords, axes)
        xyzs = add_linker_xyzs(xyzs, linker, l1_coords, axes, thetax=np.pi, trnsz=cube_l)
        xyzs = add_linker_xyzs(xyzs, linker, l1_coords, axes, thetaz=np.pi/2.0, theray=np.pi/2.0, thetax=np.pi/4.0,
                               trnsz=cube_l/2.0, trnsx=-cube_l/2.0)
        xyzs = add_linker_xyzs(xyzs, linker, l1_coords, axes, thetaz=np.pi/2.0, theray=-np.pi/2.0, thetax=-np.pi/4.0,
                               trnsz=cube_l/2.0, trnsx=cube_l/2.0)
        xyzs = add_linker_xyzs(xyzs, linker, l1_coords, axes, thetaz=np.pi/2.0, theray=-np.pi/2.0, thetaz2=-np.pi/2.0,
                               thetax=np.pi/4.0, trnsz=cube_l/2.0, trnsy=-cube_l/1.5)
        xyzs = add_linker_xyzs(xyzs, linker, l1_coords, axes, thetaz=np.pi/2.0, theray=-np.pi/2.0, thetaz2=np.pi/2.0,
                               thetax=np.pi/4.0, trnsz=cube_l/2.0, trnsy=cube_l/1.5)
        xyzs = add_metal_xyzs(xyzs, cage, linker, l1_coords, l2_coords)

        if has_ok_m_x_distances(xyzs, linker, max_mx_dist=cube_l/1.5):
            from .input_output import xyzs2xyzfile
            xyzs2xyzfile(xyzs, basename='tmp_cage')

        if is_geom_reasonable(xyzs, supress_print=True) and has_ok_m_x_distances(xyzs, linker, max_mx_dist=cube_l/1.5):
            from .input_output import xyzs2xyzfile
            xyzs2xyzfile(xyzs, basename='tmp2_cage')
            break

    if not is_geom_reasonable(xyzs):
        if not Config.suppress_print:
            print("WARNING:  It doesn't look like a sensible cage can be built")
        return None

    xyzs = uff_cage_opt(xyzs, cage, linker, n_ligands, n_metals, tmp_m_atom)

    return xyzs
