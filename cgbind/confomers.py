from rdkit import Chem
from cgbind.log import logger
import os


def extract_xyzs_from_rdkit_mol_object(mol_obj, conf_ids):
    """
    Generate xyz lists for all the conformers in conf_ids
    :param mol_obj: Molecule object
    :param conf_ids: (list) list of conformer ids to convert to xyz
    :return: (list) of xyz lists
    """
    xyzs = []

    for i in range(len(conf_ids)):
        mol_block_lines = Chem.MolToMolBlock(mol_obj, confId=conf_ids[i]).split('\n')
        mol_file_xyzs = []

        for line in mol_block_lines:
            split_line = line.split()
            if len(split_line) == 16:
                atom_label, x, y, z = split_line[3], split_line[0], split_line[1], split_line[2]
                mol_file_xyzs.append([atom_label, float(x), float(y), float(z)])

        xyzs.append(mol_file_xyzs)

    if len(xyzs) == 0:
        logger.critical('Length of conformer xyz list was 0')
        exit()

    return xyzs
