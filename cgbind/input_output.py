import os
from datetime import date
from cgbind.log import logger
from cgbind.exceptions import FileMalformatted


def xyzs2xyzfile(xyzs, filename=None, basename=None, title_line=''):
    """
    For a list of xyzs in the form e.g [[C, 0.0, 0.0, 0.0], ...] convert create a standard .xyz file

    :param xyzs: (list(list))
    :param filename: (str) Name of the generated xyz file
    :param basename: (str) Name of the generated xyz file without the file extension
    :param title_line: (str) String to print on the title line of an xyz file
    :return: (str) xyz filename
    """

    if basename is not None:
        filename = basename + '.xyz'

    if filename is None:
        logger.error('Could not print an .xyz. Filename was None')
        return None

    if xyzs is None:
        logger.error('No xyzs to print')
        return None

    if not filename.endswith('.xyz'):
        logger.error('Filename does not end with .xyz, adding')
        filename += '.xyz'

    if filename is None:
        logger.error('Filename was None')
        return filename

    with open(filename, 'w') as xyz_file:
        print(len(xyzs), '\n', title_line, sep='', file=xyz_file)
        [print('{:<3}{:^10.5f}{:^10.5f}{:^10.5f}'.format(*line), file=xyz_file) for line in xyzs]

    return filename


def xyzfile2xyzs(filename):
    """
    Convert a standard xyz file into a list of xyzs

    :param filename: (str)
    :return: (list(list))
    """
    logger.info(f'Converting {filename} to list of xyzs')

    xyzs = []

    if os.path.exists(filename) and filename.endswith('.xyz'):
        with open(filename, 'r') as xyz_file:
            xyz_lines = xyz_file.readlines()[2:]
            for line in xyz_lines:
                atom_label, x, y, z = line.split()
                xyzs.append([atom_label, float(x), float(y), float(z)])

    else:
        logger.error('Could not read .xyz file')
        return None

    if len(xyzs) == 0:
        logger.error(f'Could not read xyz lines in {filename}')
        raise FileMalformatted

    return xyzs


def molfile2xyzs(filename):
    """
    Convert a .mol file to a list of xyzs

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

    :param filename: (str)
    :return: (list(list))
    """
    xyzs = []

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
            xyzs.append([atom_label, float(x), float(y), float(z)])

    if len(xyzs) == 0:
        logger.error(f'Could not read xyz lines in {filename}')
        raise FileMalformatted

    return xyzs


def mol2file2xyzs(filename):
    """
    Convert a .mol file into a standard set of xyzs in the form e.g [[C, 0.0, 0.0, 0.0], ...]

    :param filename: (str) name of the .mol file
    :return: (lis(list)) xyzs
    """
    logger.info('Converting .mol2 file to xyzs')

    if not (os.path.exists(filename) and filename.endswith('.mol2')):
        logger.error('Could not read .mol2 file')
        raise FileMalformatted

    mol_file_lines = open(filename, 'r').readlines()

    # Get the unformatted xyzs from the .mol2 file. The atom labels will not be standard
    unformat_xyzs, xyz_block = [], False
    for n_line, line in enumerate(mol_file_lines):

        if '@' in line and xyz_block:
            break

        if xyz_block:
            try:
                atom_label, x, y, z = line.split()[1:5]
                try:
                    unformat_xyzs.append([atom_label, float(x), float(y), float(z)])
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

    xyzs = []
    for xyz_line in unformat_xyzs:
        atom_label = xyz_line[0]

        if len(atom_label) == 1:
            xyzs.append([atom_label] + xyz_line[1:])

        # e.g. Pd1 or C58
        elif atom_label[0].isalpha() and not atom_label[1].isalpha():
            xyzs.append([atom_label[0]] + xyz_line[1:])

        # e.g. Pd10
        elif atom_label[0].isalpha() and atom_label[1].isalpha():
            xyzs.append([atom_label[:2]] + xyz_line[1:])

        else:
            logger.error('Unrecognised atom type')
            return

    return xyzs
