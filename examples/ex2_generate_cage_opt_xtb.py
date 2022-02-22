# Check that autodE is available. Required to run electronic structure theory
# (EST) calculations
try:
    import autode
except ModuleNotFoundError:
    exit('ERROR: autodE (https://github.com/duartegroup/autode) must be '
         'installed to run XTB calculations. Exiting')

# Import Cage and Linker objects as well as an autode EST object for xtb
from cgbind import Linker, Cage, xtb

# Generates the same linker as in ex0
linker = Linker(name='linker_opt',
                smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                arch_name='m2l4')

# From the optimsied linker build the M2L4 cage with Pd(II) ions
cage = Cage(linker, metal='Pd', metal_charge='2')

# Optimise the cage with XTB using 2 cores. The full cage can also be optimised
# using ORCA (with DFT) in a similar way but will take several hours
cage.optimise(method=xtb, n_cores=2)

cage.print_to_file()

"""
To optimise either the linker or the full cage at the DFT level see ex3
"""
e