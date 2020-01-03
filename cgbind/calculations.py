from autode.calculation import Calculation
from autode.methods import XTB, ORCA
from cgbind.defaults import *
from cgbind.log import logger
import os
import tempfile
import shutil


def optimise(molecule, method, keywords, n_cores=1, max_core_mb=1000, cartesian_constraints=None):
    """
    Optimise a molecule

    :param molecule: (object)
    :param method: (autode.ElectronicStructureMethod)
    :param keywords: (list(str)) Keywords to use for the electronic structure calculation e.g. ['Opt', 'PBE', 'def2-SVP']
    :param n_cores: (int) Number of cores to use
    :param max_core_mb: (float)
    :param cartesian_constraints: (list(int)) List of atom ids to constrain
    :return:
    """
    logger.info('Running an optimisation calculation')

    if keywords is None:
        if method == ORCA:
            logger.warning('No keywords were set for the optimisation but an ORCA calculation was requested.\n'
                           'Using defaults.orca_low_opt_keywords..')
            keywords = orca_low_opt_keywords

        elif method == XTB:
            # No keywords are required for XTB
            pass

        else:
            logger.critical('No keywords were set for the optimisation calculation')
            exit()

    opt = Calculation(name=molecule.name + '_opt', molecule=molecule, method=method, keywords=keywords,
                      n_cores=n_cores, max_core_mb=max_core_mb, cartesian_constraints=cartesian_constraints, opt=True)
    opt.run()
    molecule.energy = opt.get_energy()
    molecule.xyzs = opt.get_final_xyzs()

    return None


def singlepoint(molecule, method, keywords, n_cores=1, max_core_mb=1000):
    """
    Run a single point energy evaluation on a molecule

    :param molecule: (object)
    :param method: (autode.ElectronicStructureMethod)
    :param keywords: (list(str)) Keywords to use for the electronic structure calculation e.g. ['Opt', 'PBE', 'def2-SVP']
    :param n_cores: (int) Number of cores to use
    :param max_core_mb: (float)
    :return:
    """
    logger.info('Running single point calculation')

    if keywords is None:
        if method == ORCA:
            logger.warning('No keywords were set for the single point but an ORCA calculation was requested.\n'
                           'Using defaults.orca_sp_keywords..')
            keywords = orca_sp_keywords

        elif method == XTB:
            # No keywords are required for XTB
            pass

        else:
            logger.critical('No keywords were set for the single-point calculation')
            exit()

    sp = Calculation(name=molecule.name + '_sp', molecule=molecule, method=method, keywords=keywords,
                     n_cores=n_cores, max_core_mb=max_core_mb)
    sp.run()
    molecule.energy = sp.get_energy()

    return None


def get_charges(molecule):
    """
    Get the partial atomic charges with XTB (tested with v. 6.2) will generate then trash a temporary directory

    :return:
    """
    logger.info('Getting charges')

    if not XTB.available:
        logger.error('Could not calculate without an XTB install')
        return None

    # Make a temporary directory
    here = os.getcwd()
    tmp_dirpath = tempfile.mkdtemp()
    os.chdir(tmp_dirpath)

    # Run the calculation
    xtb_sp = Calculation(name=molecule.name + '_xtb_sp', molecule=molecule, method=XTB, n_cores=1)
    xtb_sp.run()

    if not os.path.exists('charges'):
        logger.error('Could not get the charges from the XTB file')
        return []

    # Charges file from XTB is one value per line
    charges = [float(line.split()[0]) for line in open('charges', 'r').readlines()]

    # Remove the temporary directory
    os.chdir(here)
    shutil.rmtree(tmp_dirpath)

    if len(charges) == molecule.n_atoms:
        return charges

    else:
        logger.error('XTB failed to generate charges')
        return None
