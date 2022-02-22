import types
import functools
from subprocess import Popen, DEVNULL, PIPE, STDOUT
import os
from functools import wraps
from tempfile import mkdtemp
import shutil




def copy_func(f):
    """Based on http://stackoverflow.com/a/6528148/190597 (Glenn Maynard)"""
    g = types.FunctionType(f.__code__, f.__globals__, name=f.__name__,
                           argdefs=f.__defaults__,
                           closure=f.__closure__)
    g = functools.update_wrapper(g, f)
    g.__kwdefaults__ = f.__kwdefaults__
    return g


def work_in_tmp_dir():
    """Execute a function in a temporary directory.
    """
    def func_decorator(func):

        @wraps(func)
        def wrapped_function(*args, **kwargs):
            here = os.getcwd()

            base_dir =  None

            if base_dir is not None:
                assert os.path.exists(base_dir)

            tmpdir_path = mkdtemp(dir=base_dir)


            # Move directories and execute
            os.chdir(tmpdir_path)

            result = func(*args, **kwargs)

            os.chdir(here)

            shutil.rmtree(tmpdir_path)
            return result

        return wrapped_function
    return func_decorator




def boltzmann_distribution(energies):
    K = 300  # in kelvin
    R = 8.3144598  # gas constant J/K/mol
    RK = R * K

    normalized_energies = (2625.50*1000)*(energies - np.min(energies))
    exp = np.exp(-(normalized_energies/RK))
    weights = exp/np.sum(exp)
    return weights
