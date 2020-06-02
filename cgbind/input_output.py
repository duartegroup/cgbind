import os
from cgbind.atoms import Atom
from cgbind.log import logger
from cgbind.exceptions import FileMalformatted, CgbindCritical


def xyz_file_to_atoms(filename):
    """/
    From an .xyz file get a list of atoms

    :param filename: (str) .xyz filename
    :return: (list(Atom))
    """
    logger.info(f'Getting atoms from {filename}')

    atoms = []

    if not filename.endswith('.xyz'):
        raise FileMalformatted

    # Open the file that exists and should(!) be in the correct format
    with open(filename, 'r') as xyz_file:

        try:
            n_atoms = int(xyz_file.readline().split()[0])       # First item in an xyz file is the number of atoms

        except IndexError:
            raise FileMalformatted

        xyz_lines = xyz_file.readlines()[1:n_atoms+1]       # XYZ lines should be the following 2 + n_atoms lines

        for line in xyz_lines:

            try:
                atom_label, x, y, z = line.split()[:4]
                atoms.append(Atom(atomic_symbol=atom_label, x=float(x), y=float(y), z=float(z)))

            except (IndexError, TypeError, ValueError):
                raise FileMalformatted

    return atoms


def atoms_to_xyz_file(atoms, filename, title_line=''):
    """
    Print a standard .xyz file from a set of atoms

    :param atoms: (list(Atom))
    :param filename: (str)
    :param title_line: (str)
    """

    with open(filename, 'w') as xyz_file:
        print(len(atoms), title_line, sep='\n', file=xyz_file)
        for atom in atoms:
            x, y, z = atom.coord
            print(f'{atom.label:<3} {x:^10.5f} {y:^10.5f} {z:^10.5f}', file=xyz_file)

    return None


def get_atoms_from_file(filename):
    """Get a list of atoms from a structure file"""

    if filename.endswith('.xyz'):
        return xyzfile_to_atoms(filename)

    elif filename.endswith('.mol'):
        return molfile_to_atoms(filename)

    elif filename.endswith('.mol2'):
        return mol2file_to_atoms(filename)

    else:
        raise CgbindCritical(message='Unsupported file format')


def xyzfile_to_atoms(filename):
    """
    Convert a standard xyz file into a list of atoms

    :param filename: (str)
    :return: (list(cgbind.atoms.Atom))
    """
    logger.info(f'Converting {filename} to list of atoms')

    if not (os.path.exists(filename) and filename.endswith('.xyz')):
        logger.error('Could not read .xyz file')
        raise FileMalformatted

    atoms = []

    with open(filename, 'r') as xyz_file:
        xyz_lines = xyz_file.readlines()[2:]
        for line in xyz_lines:
            atom_label, x, y, z = line.split()
            atoms.append(Atom(atom_label, float(x), float(y), float(z)))

    if len(atoms) == 0:
        logger.error(f'Could not read xyz lines in {filename}')
        raise FileMalformatted

    return atoms


def molfile_to_atoms(filename):
    """
    Convert a .mol file to a list of atoms

    :param filename: (str)
    :return: (list(Atom))
    """
    """
    e.g. for methane:
    _____________________

     OpenBabel03272015013D

      5  4  0  0  0  0  0  0  0  0999 V2000
       -0.2783    0.0756    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        0.7917    0.0756    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
       -0.6349   -0.9294   -0.0876 H   0  0  0  0  0  0  0  0  0  0  0  0
       -0.6349    0.6539   -0.8266 H   0  0  0  0  0  0  0  0  0  0  0  0
       -0.6349    0.5022    0.9141 H   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0  0  0  0
      1  3  1  0  0  0  0
      1  4  1  0  0  0  0
      1  5  1  0  0  0  0
    M  END
    _____________________
    """

    atoms = []

    if not (os.path.exists(filename) and filename.endswith('.mol')):
        logger.error('Could not read .mol file')
        raise FileMalformatted

    with open(filename, 'r') as mol_file:
        mol_lines = mol_file.readlines()[3:]
        try:
            n_atoms = int(mol_lines[0].split()[0])

        except ValueError:
            raise FileMalformatted

        for line in mol_lines[1:n_atoms+1]:
            x, y, z, atom_label = line.split()[:4]
            atoms.append(Atom(atom_label, float(x), float(y), float(z)))

    if len(atoms) == 0:
        logger.error(f'Could not read xyz lines in {filename}')
        raise FileMalformatted

    return atoms


def mol2file_to_atoms(filename):
    """
    Convert a .mol file into a standard set of atoms

    :param filename: (str)
    :return: (lis(Atom))
    """
    logger.info('Converting .mol2 file to atoms')

    if not (os.path.exists(filename) and filename.endswith('.mol2')):
        logger.error('Could not read .mol2 file')
        raise FileMalformatted

    mol_file_lines = open(filename, 'r').readlines()

    # Get the unformatted atoms from the .mol2 file. The atom labels will not be standard
    atoms, xyz_block = [], False
    for n_line, line in enumerate(mol_file_lines):

        if '@' in line and xyz_block:
            break

        if xyz_block:
            try:
                atom_label, x, y, z = line.split()[1:5]
                try:
                    atoms.append(Atom(atom_label, float(x), float(y), float(z)))
                except TypeError:
                    logger.error('There was a problem with the .mol2 file')
                    raise FileMalformatted

            except IndexError:
                logger.error('There was a problem with the .mol2 file')
                raise FileMalformatted

        # e.g.   @<TRIPOS>ATOM
        #        1 Pd1     -2.1334  12.0093  11.5778   Pd        1 RES1   2.0000
        if '@' in line and 'ATOM' in line and len(mol_file_lines[n_line+1].split()) == 9:
            xyz_block = True

    # Fix any atom labels
    for atom in atoms:

        if len(atom.label) == 1:
            continue

        # e.g. P1 or C58
        elif atom.label[0].isalpha() and not atom.label[1].isalpha():
            atom.label = atom.label[0]

        # e.g. Pd10
        elif atom.label[0].isalpha() and atom.label[1].isalpha():
            atom.label = atom.label[:2].title()

        else:
            logger.error('Unrecognised atom type')
            raise FileMalformatted

    return atoms
