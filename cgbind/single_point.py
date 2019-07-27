from .log import logger
from .input_output import print_output
from .config import Config
from .ORCAio import get_single_point_energy


def singlepoint(mol, n_cores):

    if mol.xyzs is None:
        logger.error('Have no xyzs. Cannot run a singlepoint calculation')
        return None

    if Config.code != 'orca':
        logger.error('Only ORCA single point calculations are implemented')
        return None

    print_output('Single point calculation of', mol.name, 'Running')

    if Config.path_to_orca is None:
        logger.error('path_to_orca needs to be set for an ORCA calculation. Skipping')
        print_output('', '', 'Failed')
        return None

    energy = get_single_point_energy(mol.xyzs, mol.name, Config.sp_keywords, Config.sp_solvent, mol.charge,
                                     n_cores=n_cores)
    print_output('', '', 'Done')

    return energy
