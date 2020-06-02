import numpy as np
import itertools
import networkx as nx
from copy import deepcopy
from cgbind.exceptions import CgbindCritical
from cgbind.log import logger
from cgbind.atoms import get_metal_favoured_heteroatoms


def get_cost_metal_x_atom_interaction(x_motifs, linker, metal):
    """
    Calculate a cost based on the favourability of the M--X interaction

    :param x_motifs: ((list(Xmotif))
    :param linker: (Linker)
    :param metal: (str)
    :return:
    """

    if metal is None:
        logger.warning('Could not sort x motifs list. Metal was not specified')
        return 0

    fav_x_atoms = get_metal_favoured_heteroatoms(metal=metal)
    cost = 0

    for x_motif in x_motifs:
        for atom_id in x_motif.atom_ids:
            atom_label = linker.atoms[atom_id].label
            if atom_label in fav_x_atoms:
                cost += 10 * fav_x_atoms.index(atom_label)

    return cost


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

    return np.array(shifted_coords)


def check_x_motifs(linker=None, linker_template=None):
    if linker is None and linker_template is not None:
        if not all([motif.n_atoms == linker_template.x_motifs[0].n_atoms for motif in linker_template.x_motifs]):
            logger.critical('Found x motifs in the structure that have different number of atoms')
            exit()
        else:
            return None

    if not all([motif.n_atoms == linker_template.x_motifs[0].n_atoms for motif in linker.x_motifs]):
        logger.warning('Found x motifs in the structure that have different number of atoms')
        logger.info('Stripping the motifs with the wrong number of atoms')

        linker.x_motifs = [motif for motif in linker.x_motifs if motif.n_atoms == linker_template.x_motifs[0].n_atoms]
        logger.info(f'Now have {len(linker.x_motifs)} motifs in the linker')

    if len(linker.x_motifs) == 0:
        raise CgbindCritical(message='Have 0 Xmotifs – cannot build a cage. Is the template correct?')

    if len(linker.x_motifs) > 0:
        logger.info(f'Number of atoms in the x motifs is {linker.x_motifs[0].n_atoms}')
    return None


def powerset(s):
    """[0, 1, 2] -> [(0, 1), (1, 2) (0, 2), (0, 1 2)]"""
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(2, len(s)+1))


def is_fully_connected(atom_indexes, bonds):
    """
    Determine if a set of atoms is fully connected given a list of bonds

    :param atom_indexes: (list(int))
    :param bonds: (list(tuple))
    :return: (bool)
    """
    graph = nx.Graph()

    for atom_index in atom_indexes:
        graph.add_node(atom_index)

    for (i, j) in bonds:
        if i in atom_indexes and j in atom_indexes:
            graph.add_edge(i, j)

    return nx.is_connected(graph)


def find_x_motifs(linker):
    """
    Find the X motifs in a structure which correspond to the X atoms and their
    nearest neighbours. These may be joined if they are bonded

    :return: (list(list(int)))
    """

    def centroid_atom_distance(atom_i):
        return np.linalg.norm(linker.atoms[atom_i].coord - linker.com)

    x_motifs = []

    for donor_atom in linker.x_atoms:
        x_motif = []

        # Add all the atoms that are connected to the donor atom
        for (i, j) in linker.bonds:
            if donor_atom == i and j not in x_motif:
                x_motif.append(j)
            if donor_atom == j not in x_motif:
                x_motif.append(i)

        logger.info(f'X atom {donor_atom} had {len(x_motif)} nearest neighbours')
        x_motif.append(donor_atom)
        x_motifs.append(x_motif)

    # Get all the combinations of x motifs with length > 2 up to the total
    # number of x_motifs
    x_motif_combinations = powerset(s=deepcopy(x_motifs))

    logger.info(f'Have {len(list(powerset(s=deepcopy(x_motifs))))} groups of X'
                f' motifs to determine if they are bonded')
    for i, x_motif_group in enumerate(x_motif_combinations):
        logger.info(f'Determining if all {len(x_motif_group)} x motifs in this'
                    f' group are bonded')

        x_motif_group_atom_indexes = []
        for x_motif in x_motif_group:
            x_motif_group_atom_indexes += list(x_motif)

        if is_fully_connected(x_motif_group_atom_indexes, bonds=linker.bonds):
            logger.info(f'X-motifs are bonded')
            x_motifs.append(list(set(x_motif_group_atom_indexes)))

    logger.info(f'Found {len(x_motifs)} X motifs in the linker, '
                f'with {set([len(x) for x in x_motifs])} atoms')

    # Order the x_motifs according to the centroid – coord
    # distance: smallest -> largest
    sorted_x_motifs_ids = [sorted(list(x_motif), key=centroid_atom_distance)
                           for x_motif in x_motifs]

    return [Xmotif(atom_ids=motif, coords=[linker.atoms[i].coord for i in motif]) for motif in sorted_x_motifs_ids]


def get_maximally_connected_x_motifs(x_motifs, x_atoms):
    """
    Given a list of Xmotifs find those that are maximally connected, i.e. the ones that contain all the donor atoms
    but also are the largest in size

    :param x_motifs: (list(cgbind.x_motifs.Xmotif)
    :param x_atoms: (list(int))
    :return:
    """

    #                      X motif lengths sorted from high to low
    for x_motif_length in reversed(sorted(set([len(x) for x in x_motifs]))):

        new_x_motifs = [x for x in x_motifs if len(x) == x_motif_length]

        # Add all the atom ids of the xmotifs to a single list
        x_motifs_atoms = []
        for x_motif in new_x_motifs:
            x_motifs_atoms += x_motif.atom_ids

        # All the donor (X) atoms need to be in the full list
        if all(x_atom in x_motifs_atoms for x_atom in x_atoms):
            logger.info(f'Returning {len(new_x_motifs)} Xmotifs each with {len(new_x_motifs[0])} atoms')
            return new_x_motifs

    logger.critical('Could not find a set of x motifs of the same length with'
                    ' all the donor atoms')
    raise CgbindCritical


class Xmotif:

    def __len__(self):
        return len(self.atom_ids)

    def __init__(self, atom_ids, coords):
        """
        Xmotif. e.g. CNC in a classic M2L4 pyridyl based linker

        :param atom_ids: (list(int))
        :param coords: (list(np.ndarray))
        """

        self.atom_ids = atom_ids               #: (list(int)) Atom ids of the X motif
        self.coords = coords                   #: (list(np.ndarray)) Coordinates of the x motif atoms
        self.n_atoms = len(coords)             #: (int) Number of atoms in the X motif
        self.shift_vec = None                  #: (np.ndarray) Vector to shift by when expanding/contracting
        self.r = None                          #: (float) Current distance from the cage centroid to the Xmotif centroid
        self.norm_shift_vec = None             #: (np.ndarray) Normalised shift vector
