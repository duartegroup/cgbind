from cgbind.defaults import *
from cgbind.log import logger
from cgbind.exceptions import RequiresAutodE


def optimise(molecule, method, keywords, n_cores=1, cartesian_constraints=None):
    """
    Optimise a molecule

    :param molecule: (object)
    :param method: (autode.ElectronicStructureMethod)
    :param keywords: (list(str)) Keywords to use for the electronic structure calculation e.g. ['Opt', 'PBE', 'def2-SVP']
    :param n_cores: (int) Number of cores to use
    :param cartesian_constraints: (list(int)) List of atom ids to constrain
    :return:
    """
    logger.info('Running an optimisation calculation')

    try:
        from autode.calculation import Calculation
        from autode.wrappers.XTB import xtb
        from autode.wrappers.ORCA import orca

    except ModuleNotFoundError:
        logger.error('autode not found. Calculations not available')
        raise RequiresAutodE

    if keywords is None:
        if method == orca:
            logger.warning('No keywords were set for the optimisation but an ORCA calculation was requested.\n'
                           'Using defaults.orca_low_opt_keywords..')
            keywords = orca_low_opt_keywords

        elif method == xtb:
            # No keywords are required for XTB
            pass

        else:
            logger.critical('No keywords were set for the optimisation calculation')
            exit()

    opt = Calculation(name=molecule.name + '_opt', molecule=molecule, method=method, keywords_list=keywords,
                      n_cores=n_cores, cartesian_constraints=cartesian_constraints, opt=True)
    opt.run()
    molecule.energy = opt.get_energy()
    molecule.set_atoms(atoms=opt.get_final_atoms())

    return None


def singlepoint(molecule, method, keywords, n_cores=1):
    """
    Run a single point energy evaluation on a molecule

    :param molecule: (object)
    :param method: (autode.ElectronicStructureMethod)
    :param keywords: (list(str)) Keywords to use for the electronic structure calculation e.g. ['Opt', 'PBE', 'def2-SVP']
    :param n_cores: (int) Number of cores to use
    :return:
    """
    logger.info('Running single point calculation')

    try:
        from autode.calculation import Calculation
        from autode.wrappers.XTB import xtb
        from autode.wrappers.ORCA import orca

    except ModuleNotFoundError:
        logger.error('autode not found. Calculations not available')
        raise RequiresAutodE

    if keywords is None:
        if method == orca:
            logger.warning('No keywords were set for the single point but an ORCA calculation was requested.\n'
                           'Using defaults.orca_sp_keywords..')
            keywords = orca_sp_keywords

        elif method == xtb:
            # No keywords are required for XTB
            pass

        else:
            logger.critical('No keywords were set for the single-point calculation')
            exit()

    sp = Calculation(name=molecule.name + '_sp', molecule=molecule, method=method, keywords_list=keywords,
                     n_cores=n_cores)
    sp.run()
    molecule.energy = sp.get_energy()

    return None


def get_charges(molecule):
    """
    Get the partial atomic charges with XTB (tested with v. 6.2) will generate then trash a temporary directory

    :return:
    """
    logger.info('Getting charges')

    try:
        from autode.calculation import Calculation
        from autode.wrappers.XTB import xtb
        from autode.exceptions import MethodUnavailable

    except ModuleNotFoundError:
        logger.error('autode not found. Calculations not available')
        raise RequiresAutodE

    # Run the calculation
    try:
        xtb_sp = Calculation(name=molecule.name + '_xtb_sp', molecule=molecule, method=xtb, n_cores=1)
        xtb_sp.run()

        charges = xtb_sp.get_atomic_charges()

    except MethodUnavailable:
        logger.error('Could not calculate without an XTB install')
        return None

    if len(charges) == molecule.n_atoms:
        return charges

    else:
        logger.error('XTB failed to generate charges')
        return None
