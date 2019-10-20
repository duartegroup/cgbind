# Import all the functions from cgbind
from cgbind import *

# Set the optimisation method to the 'min_dft' (PBE-D3BJ/def2-SVP) level available in ORCA
Config.code = 'orca'
Config.path_to_orca = '/Users/tom/opt/orca_4_1/orca'
Config.opt_method = 'min_dft'

# Generate a cage.substrate complex and optimise just the substrate inside the fixed cage geometry
gen_cage_subst_complex('L-1', 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                       'q1', 'O=C1C=CC(C=C1)=O',
                       opt_cage=False,
                       opt_linker=False,
                       opt_substrate=False,
                       opt_cage_subst=False,
                       fix_cage_geom=False,
                       arch_name='m2l4')

"""
To calculate binding affinities see ex5
"""
