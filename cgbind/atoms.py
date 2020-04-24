from cgbind.log import logger
import numpy as np


class Atom:

    def __init__(self, atomic_symbol, x, y, z):

        self.label = atomic_symbol
        self.coord = np.array([float(x), float(y), float(z)])


heteroatoms = ['O', 'N', 'S', 'P', 'F', 'Cl']

metals = ['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga',
          'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Cs', 'Ba', 'La', 'Ce',
          'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os',
          'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
          'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl',
          'Mc', 'Lv']

avg_bond_lengths = {
    'HH': 0.74,
    'CC': 1.54,
    'NN': 1.45,
    'OO': 1.48,
    'FF': 1.42,
    'ClCl': 1.99,
    'II': 2.67,
    'CN': 1.47,
    'CO': 1.43,
    'CS': 1.82,
    'CF': 1.35,
    'CCl': 1.77,
    'CBr': 1.94,
    'CI': 2.14,
    'HC': 1.09,
    'HN': 1.01,
    'HO': 0.96,
    'HBr': 1.41,
    'HCl': 1.27,
    'HI': 1.61
}

atomic_masses = {
    'H': 1.01,
    'B': 10.81,
    'C': 12.01,
    'N': 14.01,
    'O': 16.0,
    'F': 19.0,
    'Si': 28.09,
    'P': 32.07,
    'S': 32.07,
    'Cl': 35.45,
    'Br': 79.90,
    'I': 126.90
}


# Van der Waals radii in Å taken from http://www.webelements.com/periodicity/van_der_waals_radius/

vdw_radii = {
    'H': 1.20,
    'He': 1.40,
    'Li': 1.82,
    'Be': 1.53,
    'B': 1.92,
    'C': 1.70,
    'N': 1.55,
    'O': 1.52,
    'F': 1.47,
    'Ne': 1.54,
    'Na': 2.27,
    'Mg': 1.73,
    'Al': 1.84,
    'Si': 2.10,
    'P': 1.80,
    'S': 1.80,
    'Cl': 1.75,
    'Ar': 1.88,
    'K': 2.75,
    'Ca': 2.31,
    'Ni': 1.63,
    'Cu': 1.40,
    'Zn': 1.39,
    'Ga': 1.87,
    'Ge': 2.11,
    'As': 1.85,
    'Se': 1.90,
    'Br': 1.85,
    'Kr': 2.02,
    'Rb': 3.03,
    'Sr': 2.49,
    'Pd': 1.63,
    'Ag': 1.72,
    'Cd': 1.58,
    'In': 1.93,
    'Sn': 2.17,
    'Sb': 2.06,
    'Te': 2.06,
    'I': 1.98,
    'Xe': 2.16,
    'Cs': 3.43,
    'Ba': 2.49,
    'Pt': 1.75,
    'Au': 1.66}

atoms = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
         'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
         'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn' 'Sb', 'Te',
         'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm',
         'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
         'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
         'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

# TODO: Make this less empirical... e.g. tabulate all M-L bond enthalpies using CCSDT//DFT
metals_and_favored_heteroatoms = {
    'Pd': ['P', 'N', 'S', 'O', 'F', 'Cl']
}


def get_max_valency(atom):

    if atom.label in ['H', 'F']:
        return 1

    if atom.label in ['B', 'Al', 'Ga']:
        return 4

    elif atom.label in ['C', 'Si', 'Ge']:
        return 4

    elif atom.label in ['N', 'P', 'As']:
        return 4

    elif atom.label in ['O', 'S', 'Se']:
        return 3

    elif atom.label == 'F':
        return 2

    else:
        logger.warning('Could not find a maximum valency, returning 6')
        return 6


def get_metal_favoured_heteroatoms(metal):
    """
    For a metal return the list of favoured donor atoms from most -> least favoured
    :param metal: (str) Atomic symbol
    :return: (list(str))
    """
    assert type(metal) == str

    try:
        return metals_and_favored_heteroatoms[metal]

    except KeyError:
        logger.warning(f'Favoured {metal}--donor interactions are not tabulated')
        return heteroatoms


def get_atomic_number(atom):
    """
    Get the atomic number of an atom given by its atomic symbol e.g. Kr --> 36

    :param atom: (cgbind.atoms.Atom)
    :return: (int) Atomic number
    """

    try:
        return atoms.index(atom.label) + 1

    except ValueError:
        logger.warning('Could not get atomic number. Returning 6..')
        return 6


def get_atomic_mass(atom):
    """
    Get the atomic mass in amu given the atomic symbol e.g. C --> 12.01

    :param atom: (cgbind.atoms.Atom)
    :return: (int) Atomic number
    """

    if atom.label in atomic_masses.keys():
        atom_mass = atomic_masses[atom.label]
    else:
        logger.warning(f"Couldn't find the atomic mass for {atom.label}. Guessing at 10")
        atom_mass = 10

    return atom_mass


def get_vdw_radii(atom):
    """
    Get the van der Wall radii of an atom givne its atomic symbol e.g. Kr --> ~2 Å

    :param atom: (cgbind.atoms.Atom)
    :return: (float) VdW radius in Å
    """

    if atom.label in vdw_radii.keys():
        vdv_radii = vdw_radii[atom.label]
    else:
        logger.error(f"Couldn't find the VdV radii for {atom.label}. Guessing at 1.5")
        vdv_radii = 1.5

    return vdv_radii
