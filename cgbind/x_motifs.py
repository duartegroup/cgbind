import numpy as np
from copy import deepcopy
from cgbind.log import logger


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


def check_x_motifs(linker, linker_template):
    if not all([motif.n_atoms == linker_template.x_motifs[0].n_atoms for motif in linker.x_motifs]):
        logger.warning('Found x motifs in the structure that have different number of atoms')
        logger.info('Stripping the motifs with the wrong number of atoms')

        linker.x_motifs = [motif for motif in linker.x_motifs if motif.n_atoms == linker_template.x_motifs[0].n_atoms]
        logger.info(f'Now have {len(linker.x_motifs)} motifs in the linker')

    if len(linker.x_motifs) > 0:
        logger.info(f'Number of atoms in the x motifs is {linker.x_motifs[0].n_atoms}')
    return None


def find_x_motifs(linker):
    """
    Find the X motifs in a structure which correspond to the X atoms and their nearest neighbours. These motifs will
    be searched in a linker object and the RMSD minimised

    :return: (list(list(int)))
    """

    def centroid_atom_distance(atom_i):
        return np.linalg.norm(linker.coords[atom_i] - linker.centroid)

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

    # Combine x_motifs that are bonded, thus should be considered a single motif
    bonded_x_motif_sets = []
    for n, x_motif_i in enumerate(x_motifs):
        bonded_x_motif = x_motif_i.copy()

        # Loop over all other motifs that are not the same
        for m, x_motifs_j in enumerate(x_motifs):
            if n != m:

                for (i, j) in linker.bonds:
                    # Check that there is a bond between x_motif i and x_motif j
                    if (i in bonded_x_motif and j in x_motifs_j) or (j in bonded_x_motif and i in x_motifs_j):
                        bonded_x_motif += x_motifs_j
                        break

        bonded_x_motif_sets.append(set(bonded_x_motif))

    # Add the largest set to the bonded_x_motifs as a list. Some motifs will be missed due to the oder in which
    # they're added

    largest_unique_bonded_x_motif_sets = []
    # Sort the sets according to size, then don't add identical sets or subsets
    for x_motif in reversed(sorted(bonded_x_motif_sets, key=len)):

        unique = True
        for unique_x_motif in largest_unique_bonded_x_motif_sets:

            # If the motif is already in the unique list then don't append, nor if the motif is a subset
            if x_motif == unique_x_motif or x_motif.issubset(unique_x_motif):
                unique = False
                break

        if unique:
            largest_unique_bonded_x_motif_sets.append(x_motif)

    logger.info(f'Found {len(largest_unique_bonded_x_motif_sets)} X motifs in the linker')

    # Order the x_motifs according to the centroid – coord distance: smallest -> largest
    sorted_x_motifs_ids = [sorted(list(x_motif), key=centroid_atom_distance)
                           for x_motif in largest_unique_bonded_x_motif_sets]

    return [Xmotif(atom_ids=motif, coords=[linker.coords[i] for i in motif]) for motif in sorted_x_motifs_ids]


class Xmotif:

    def __init__(self, atom_ids, coords):

        self.atom_ids = atom_ids
        self.coords = coords
        self.n_atoms = len(coords)
        self.shift_vec = None
        self.r = None
        self.norm_shift_vec = None
