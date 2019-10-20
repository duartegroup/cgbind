# Import all the functions from cgbind
from cgbind import *

# Set the optimisation method to the min_dft : PBE-D3BJ/def2-SVP level available in ORCA
Config.code = 'orca'
Config.opt_method = 'min_dft'

# Set the number of cores to 4
Config.n_cores = 4

# Note: This will work if the ORCA directory is in your $PATH, otherwise set
# Config.path_to_orca = '/path/to/orca/install/orca'

# Generate a cage from a linker whose geometry is optimised at the PBE-D3BJ/def2-SVP level. This should take ~10 mins
# to complete
gen_cage(linker_name='L-1',
         linker_smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
         opt_linker=True,
         opt_cage=False)

"""
To generate multiple cages in parallel see ex3
"""
