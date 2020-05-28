from cgbind import geom
from cgbind.config import Config
from cgbind.constants import Constants
from cgbind.linker import Linker
from cgbind.cage import Cage
from cgbind.substrate import Substrate
from cgbind.cage_subt import CageSubstrateComplex
from cgbind.templates import Template

__version__ = '1.0.0'

__all__ = ['Config',
           'Constants',
           'geom',
           'Linker',
           'Cage',
           'Substrate',
           'CageSubstrateComplex']

try:
    from autode.wrappers.MOPAC import mopac
    from autode.wrappers.ORCA import orca
    from autode.wrappers.XTB import xtb
    __all__ += ['mopac', 'orca', 'xtb']

except ModuleNotFoundError:
    pass
