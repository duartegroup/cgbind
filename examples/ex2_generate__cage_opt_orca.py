# Import all the functions from cgbind
from cgbind import *

# Set the optimisation method to the min_dft : PBE-D3BJ/def2-SVP level available in ORCA
Config.code = 'orca'
Config.path_to_orca = '/Users/tom/opt/orca_4_1/orca'
Config.opt_method = 'min_dft'

# Generate a cage from a linker whose geometry is optimised at the PBE-D3BJ/def2-SVP level, then optimise the full cage
gen_cage(linker_name='L-1',
         linker_smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
         opt_linker=True,
         opt_cage=True)

"""
To generate multiple cages in parallel see ex3
"""
