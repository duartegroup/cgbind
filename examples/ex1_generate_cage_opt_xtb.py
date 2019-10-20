# Import all the functions from cgbind
from cgbind import *

# Set the optimisation method to the tight binding GFN-XTB hamiltonian from XTB
Config.code = 'xtb'


# Note: This will work if the xtb/bin directory is in your $PATH, otherwise set
# Config.path_to_xtb = '/path/to/xtb/install/bin/xtb'

# Generate a cage from a linker whose geometry is optimised at the XTB level of theory
gen_cage(linker_name='L-1',
         linker_smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
         opt_linker=True,
         opt_cage=False,
         arch_name='m2l4')

"""
To optimise either the linker or the full cage at the DFT level see ex2
"""
