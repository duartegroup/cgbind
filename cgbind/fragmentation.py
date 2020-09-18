from cgbind.geom import get_rot_mat_kabsch
from cgbind.geom import rotation_matrix
from cgbind.bonds import get_avg_bond_length
from scipy.spatial import distance_matrix
from scipy.optimize import minimize, Bounds
from copy import deepcopy
import networkx as nx
import numpy as np


def get_frag_minimised_linker(linker, fragments, template_linker, x_motifs):
    """
    From a linker with a set of fragments, fit the fragments containing
    x-motifs that is minimused with respect to dr

    :param linker:
    :param fragments:
    :param template_linker:
    :param x_motifs:
    :return:
    """
    if len(fragments) != len(x_motifs) + 1:
        # Stitch the non x_motif containing fragments back together
        raise NotImplementedError

    # Minimise the cost function as a function of dr in Å
    opt = minimize(frag_linker,
                   x0=np.array([0.0]),
                   args=(linker, fragments, template_linker, x_motifs, True),
                   method='L-BFGS-B', tol=1e-3,
                   bounds=Bounds(-100.0, 10.0))

    # Create a linker conformer from the minimisation
    conformer = deepcopy(linker)
    conformer.dr = opt.x[0]
    conformer.x_motifs = x_motifs

    coords = frag_linker(conformer.dr, linker, fragments, template_linker,
                         x_motifs, return_coords=True)
    conformer.set_atoms(coords=coords)
    conformer.cost = opt.fun
    # linker.print_xyz_file()

    return conformer


def frag_linker(dr, linker, fragments, template_linker, x_motifs,
                return_energy=False, return_coords=False):
    """
    From a linker with a set of fragments, fit the fragments containing
    x-motifs using a defined dr to the template, then minimise the remaining
    fragments under a

    V = bonds + repulsion

    potential with respect to rotation and translation of the fragment.

    :param dr: (float) Shift value for the template
    :param linker: (cgbind.linker.Linker)
    :param fragments: (list(list(int))) List of atom indexes of each fragment
    :param template_linker: (cgbind.templates.Linker)
    :param x_motifs: (list(cgbind.x_motifs.Xmotif))
    :param return_energy: (bool) Only return the energy of the linker
    :param return_coords: (bool) Return the coordinates of the linker
    :return: (float or np.ndarray or None)
    """
    coords = np.array(linker.get_coords(), copy=True)

    fitted_idxs = []
    for i, fragment in enumerate(fragments):
        if is_x_motif_and_fit(fragment, template_linker, coords, dr, x_motifs):
            fitted_idxs.append(i)

    # Fitted coordinates are now fixed in space
    fixed_idxs = [idx for i in fitted_idxs for idx in fragments[i]]

    # Deal with the remaining, not fitted fragment
    fragment = [frag for i, frag in enumerate(fragments)
                if i not in fitted_idxs][0]

    # Minimise the energy with respect to the remaining fragments
    fragment_idxs = np.array(fragment, dtype=np.int)
    ij_bonds = get_cross_bonds(linker, fragment, fixed_idxs)

    opt = minimize(rotated_fragment_linker_energy,
                   x0=np.random.uniform(size=7),
                   args=(coords, fragment_idxs, fixed_idxs, ij_bonds),
                   method='BFGS')

    # print(opt)
    energy = rotated_fragment_linker_energy(opt.x, coords, fragment_idxs,
                                            fixed_idxs,
                                            ij_bonds,
                                            modify=True)
    if return_energy:
        return energy

    if return_coords:
        return coords

    return None


def rotated_fragment_linker_energy(x, coords, idxs_to_shift, fixed_idxs,
                                   ij_bonds,
                                   modify=False):
    """
    Rotate a fragment of the linker using a numpy array defining the shift
    vector (length 3) and the rotation (length 4)

    :param x:
    :param coords:
    :param idxs_to_shift:
    :param fixed_idxs:
    :param ij_bonds:
    :param modify: (bool) Modify the coordinates in place
    :return:
    """
    if not modify:
        coords = np.array(coords, copy=True)

    coords[idxs_to_shift] += x[:3]

    rot_mat = rotation_matrix(axis=x[3:6],
                              theta=x[6])

    coords[idxs_to_shift] = rot_mat.dot(coords[idxs_to_shift].T).T

    energy = linker_energy(coords_i=coords[fixed_idxs],
                           coords_j=coords[idxs_to_shift],
                           ij_bonds=ij_bonds)
    return energy


def linker_energy(coords_i, coords_j, ij_bonds):
    """
    Calculate the energy(/cost) of this current fragmented linker structure
    by calculating the repulsive energy as

    V = sum_ij 1/roj^2 + sum_nm 10 * (d_nm - d_nm^ideal)^2

    where i and j index atoms from coords i and j respectively and n and m
    index bonds that are present in the non fragmented linker

    :param coords_i:
    :param coords_j:
    :param ij_bonds:
    :return:
    """

    dist_mat = distance_matrix(coords_i, coords_j)
    repulsion = np.sum(np.power(dist_mat, -2))

    bonds = sum(10*(dist_mat[i, j] - dist)**2
                for (i, j), dist in ij_bonds.items())

    return repulsion + bonds


def get_cross_bonds(linker, fragment, fixed):
    """
    For a linker comprised of a fixed component and a movable fragment return
    the list of indexes for bonds between the two. e.g.

    fragment = [0, 4, 5, 9]
    fixed = [1, 6, 7, 10, 11, 12]

    where in the linker graph there is a bond between atoms 4 and 10 between
    the fragment and fixed component then this function would return
    {(1, 3): 1.5} if atoms 4 and 10 were both carbons

    :param linker:
    :param fragment: (list(int)) Atom indexes
    :param fixed: (list(int)) Atom indexes
    :return: (dict(tuple(int)))
    """

    i_idxs, j_idxs = [], []

    for (i, j) in get_cross_fragment_bonds(linker, fragment):
        if i in fixed:
            i_idxs.append(fixed.index(i))
            j_idxs.append(fragment.index(j))

        if j in fixed:
            i_idxs.append(fixed.index(j))
            j_idxs.append(fragment.index(i))

    def dist(i, j):
        return get_avg_bond_length(linker.atoms[i].label,
                                   linker.atoms[j].label)

    return {pair: dist(*pair) for pair in zip(i_idxs, j_idxs)}


def get_cross_fragment_bonds(linker, fragment):
    """
    Get a list of bonds (as a pair of atom indexes) that span between a given
    fragment and the rest of the graph

    :param linker: (cgbind.linker.Linker)
    :param fragment: (list(int))
    :return: (list(tuple(int)))
    """
    cross_bonds = []

    for (i, j) in linker.bonds:

        if i in fragment and j not in fragment:
            cross_bonds.append((i, j))

        if j in fragment and i not in fragment:
            cross_bonds.append((i, j))

    return cross_bonds


def is_x_motif_and_fit(fragment, template_linker, coords, dr, x_motifs):
    """
    If this fragment contains an x motif then fit the coordinates will
    modify coords in place

    :param fragment:
    :param template_linker:
    :param coords:
    :param dr:
    :param x_motifs:
    :return:
    """

    for j, x_motif in enumerate(x_motifs):
        if all(atom_idx in fragment for atom_idx in x_motif.atom_ids):
            # Shifted template coordinates for this linker
            t_coords = np.array(template_linker.x_motifs[j].coords,
                                copy=True)
            t_coords += dr * template_linker.x_motifs[j].norm_shift_vec

            x_coords = coords[np.array(x_motif.atom_ids, dtype=np.int)]
            # Centre on the average x, y, z coordinate
            x_centre = np.average(x_coords, axis=0)
            x_coords -= x_centre

            # Also centre the template coordinates
            t_centre = np.average(t_coords, axis=0)
            t_coords -= t_centre

            rot_mat = get_rot_mat_kabsch(p_matrix=x_coords,
                                         q_matrix=t_coords)

            # Shift to the origin, rotate the fragment atoms then shift
            # back to the correct position
            f_idxs = np.array(fragment, dtype=np.int)

            coords[f_idxs] -= x_centre
            coords[f_idxs] = rot_mat.dot(coords[f_idxs].T).T
            coords[f_idxs] += t_centre

            return True

    return False


def fragment_graph(bonds, graph, x_motifs):
    """
    Fragment a graph across single bonds that aren't within X-motifs

    :param bonds: (list(tuple(int)))
    :param graph: (nx.Graph)
    :param x_motifs: (cgbind.x_motifs.Xmotif)
    :return: (int) Number of deletions made to the graph
    """
    n_deleted = 0

    for (i, j) in bonds:

        # Already fragmented over this bond
        if (i, j) not in graph.edges:
            continue

        # Don't split over a bond that is within an x-motif
        if any((i in x_motif.atom_ids and j in x_motif.atom_ids)
               for x_motif in x_motifs):
            continue

        graph.remove_edge(i, j)
        components = list(nx.connected_components(graph))

        # Only consider single bonds that are not in rings i.e. the
        # edge deletion should generate at least another component
        if len(components) == 1:
            graph.add_edge(i, j)
            continue

        # Don't split into very small fragments
        if any(len(frag) < 6 for frag in components):
            graph.add_edge(i, j)
            continue

        n_deleted += 1

    return n_deleted
