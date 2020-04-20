import numpy as np
import itertools
from copy import deepcopy
from cgbind.exceptions import CgbindCritical
from cgbind.log import logger
from cgbind.atoms import get_metal_favoured_heteroatoms


def sort_x_motifs(x_motifs_list, linker, metal):
    """
    Sort a list of X motifs by the favourability of the M--X interaction

    :param x_motifs_list: (list((list(Xmotif)))
    :param linker: (Linker)
    :param metal: (str)
    :return:
    """
    logger.info('Sorting the X motif list by the best M--X interactions')

    if metal is None:
        logger.warning('Could not sort x motifs list. Metal was not specified')
        return x_motifs_list

    fav_x_atoms = get_metal_favoured_heteroatoms(metal=metal)
    x_motifs_list_and_cost = {}

    for x_motifs in x_motifs_list:
        cost = 0

        for x_motif in x_motifs:
            for atom_id in x_motif.atom_ids:
                atom_label = linker.atoms[atom_id].label
                if atom_label in fav_x_atoms:
                    cost += fav_x_atoms.index(atom_label)

        x_motifs_list_and_cost[x_motifs] = cost

    return sorted(x_motifs_list_and_cost, key=x_motifs_list_and_cost.get)


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


def find_x_motifs(linker):
    """
    Find the X motifs in a structure which correspond to the X atoms and their nearest neighbours. These may be joined
    if they are bonded

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

    # Get all the combinations of x motifs with length > 2 up to the total number of x_motifs
    x_motif_combinations = powerset(s=deepcopy(x_motifs))

    logger.info(f'Have {len(list(powerset(s=deepcopy(x_motifs))))} groups of X motifs to determine if they are bonded')
    for i, x_motif_group in enumerate(x_motif_combinations):
        logger.info(f'Determining if all {len(x_motif_group)} x motifs in this group are bonded')

        inter_x_motif_bonds = 0
        bonded = False

        for (x_motif_i, x_motif_j) in itertools.combinations(x_motif_group, 2):
            for atom_index_i in x_motif_i:
                for atom_index_j in x_motif_j:
                    if (atom_index_i, atom_index_j) in linker.bonds or (atom_index_j, atom_index_i) in linker.bonds:

                        # If atoms i and j are bonded then we can have this x motif
                        bonded = True
                        break

                if bonded:
                    break

            if bonded:
                inter_x_motif_bonds += 1

        # All x motifs in the group need to be bonded to each other, with the n. bonds -1 e.g. 3 Xmotifs need 2 bonds
        if inter_x_motif_bonds == len(x_motif_group) - 1:

            # Construct the set of all atoms in the bonded Xmotif (set so no repeated atoms)
            bonded_x_motif = []
            for x_motif in x_motif_group:
                bonded_x_motif += list(x_motif)

            x_motifs.append(list(set(bonded_x_motif)))

    logger.info(f'Found {len(x_motifs)} X motifs in the linker, with {set([len(x) for x in x_motifs])} atoms')

    # Order the x_motifs according to the centroid – coord distance: smallest -> largest
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

    logger.critical('Could not find a set of x motifs of the same length with all the donor atoms')
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
