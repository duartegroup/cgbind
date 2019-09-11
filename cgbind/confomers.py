from rdkit import Chem
from cgbind.log import logger
import os


def gen_conformer_mol_files(mol):
    logger.info('Generating conformer mol files for {}'.format(mol.name))

    try:
        for i in range(len(mol.conf_ids)):
            Chem.MolToMolFile(mol.mol_obj, mol.conf_filenames[i], confId=mol.conf_ids[i])
        return 0
    except AttributeError or RuntimeError:
        logger.error('Couldn\'t convert mol objects to .mol files for {}'.format(mol.name))
        return 1


def confomer_mol_files_to_xyzs(filenames, n_atoms):
    """
    For list of .mol files extract the coordinates into the standard .xyz format
    :return: List of xyzs in the form [file1_xyzs, file2_xyzs ...], where file1_xyzs = [[X, 0.0, 0.0, 0.0]...]
    """
    logger.info('Converting .mol files into xyz lists')

    xyzs = []

    for filename in filenames:
        mol_file_lines = [line for line in open(filename, 'r')]
        mol_file_xyzs = []
        for line in mol_file_lines:
            split_line = line.split()
            if len(split_line) == 16:
                atom_label, x, y, z = split_line[3], split_line[0], split_line[1], split_line[2]
                mol_file_xyzs.append([atom_label, float(x), float(y), float(z)])

        if len(mol_file_xyzs) != n_atoms:
            logger.error('There was a problem converting .mol to xyz list. n_atoms doesn\'t match')
            return
        else:
            os.remove(filename)
            xyzs.append(mol_file_xyzs)

    return xyzs
