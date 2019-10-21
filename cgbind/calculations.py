from autode.calculation import Calculation
from cgbind.log import logger


def optimise(molecule, method, keywords, n_cores=1, max_core_mb=1000, cartesian_constraints=None):
    logger.info('Running single point calculation')

    opt = Calculation(name=molecule.name + '_opt', molecule=molecule, method=method, keywords=keywords,
                      n_cores=n_cores, max_core_mb=max_core_mb, cartesian_constraints=cartesian_constraints, opt=True)
    opt.run()
    molecule.energy = opt.get_energy()
    molecule.xyzs = opt.get_final_xyzs()

    return None


def singlepoint(molecule, method, keywords, n_cores=1, max_core_mb=1000):
    logger.info('Running single point calculation')

    sp = Calculation(name=molecule.name + '_sp', molecule=molecule, method=method, keywords=keywords,
                     n_cores=n_cores, max_core_mb=max_core_mb)
    sp.run()
    molecule.energy = sp.get_energy()

    return None
