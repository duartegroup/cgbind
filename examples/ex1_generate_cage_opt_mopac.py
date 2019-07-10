# Import all the functions from cgbind
from cgbind import *

# Set the optimisation method to the semi-empirical PM7 method available in MOPAC
Config.code = 'mopac'
Config.opt_method = 'pm7'
Config.path_to_mopac = '/Users/tom/opt/mopac/MOPAC2016.exe'
Config.path_to_mopac_licence = '/Users/tom/opt/mopac/'

# Generate a cage from a linker whose geometry is optimised at the PM7 level of theory
gen_cage(linker_name='L-1',
         linker_smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
         opt_linker=True,
         opt_cage=False)

"""
To optimise either the linker or the full cage at the DFT level Config.path_to_orca needs to be set. See ex2
"""
