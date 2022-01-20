import types
import functools
import shutil
import os
from time import time
from subprocess import Popen, PIPE
from cgbind.exceptions import FileMalformatted
from cgbind.config import Config
from cgbind.log import logger
from cgbind.input_output import xyz_file_to_atoms


def copy_func(f):
    """Based on http://stackoverflow.com/a/6528148/190597 (Glenn Maynard)"""
    g = types.FunctionType(f.__code__, f.__globals__, name=f.__name__,
                           argdefs=f.__defaults__,
                           closure=f.__closure__)
    g = functools.update_wrapper(g, f)
    g.__kwdefaults__ = f.__kwdefaults__
    return g


def fast_xtb_opt(molecule, n_cycles=5, n_cores=None):
    """
    Run an optimisation on a molecule using XTB for a defined number of
    optimisation cycles

    :param molecule: (cgbind.molecules.BaseStruct)
    :param n_cycles: (int) Number of optimisation cycles to perform
    :param n_cores: (int or None)
    :return: (cgbind.molecules.BaseStruct)
    """
    start_time = time()

    if shutil.which('xtb') is None:
        logger.error('Could not optimise - no XTB install')
        return

    if molecule.n_atoms == 0:
        logger.error('Could not optimise - no atoms')
        return

    if n_cores is None:
        n_cores = Config.n_cores
    else:
        n_cores = int(n_cores)

    logger.info(f'Optimising with {n_cores} for a maximum of {n_cycles}')
    molecule.print_xyz_file(filename='mol.xyz')

    xtb_opt = Popen(['xtb', 'mol.xyz', '--opt', 'crude',
                     '--cycles', str(int(n_cycles))], stdout=PIPE, stderr=PIPE)
    _, err = xtb_opt.communicate()

    if not os.path.exists('xtbopt.xyz'):
        logger.error('XTB optimisation failed - no output generated.'
                     f' Error: {err}')
        return

    # Try to set the optimised atoms
    try:
        molecule.set_atoms(atoms=xyz_file_to_atoms('xtbopt.xyz'))

    except FileMalformatted:
        logger.error('Incorrectly formatted xyz output')

    # Clean up the generated files
    possible_files = ['mol.xyz', 'charges', 'wbo', 'xtbopt.log',
                      'xtbopt.xyz', 'xtbrestart', 'xtbtopo.mol']

    for filename in possible_files:
        if os.path.exists(filename):
            os.remove(filename)

    logger.info(f'XTB optimisation successful in {time() - start_time:.1f} s. '
                f'Optimised atoms set.')
    return None
