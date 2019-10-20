# Import all the functions from cgbind
from cgbind import *

# Set the optimisation method to the 'min_dft' (PBE-D3BJ/def2-SVP) level available in ORCA
Config.code = 'orca'
Config.path_to_orca = '/Users/tom/opt/orca_4_1/orca'
Config.opt_method = 'min_dft'

# Set the keywords to use for a single point calculation in ORCA, used to calculate binding affinities
Config.sp_keywords = ['SP', 'PBE', 'RI', 'def2-SVP']

# To accelerate this example the two calculations are performed in parallel
Config.n_cores = 2

# Calculate binding affinities for cages generated from L-1 and L-2
calc_binding_affinities({'L-1': 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                         'L-2': 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=N2)=CC=CN=C1'},
                        {'quinone': 'O=C1C=CC(C=C1)=O'},
                        opt_linker=False,
                        opt_cage=True,
                        sp_cage=True,
                        opt_substrate=True,
                        sp_substrate=True,
                        opt_cage_subst=True,
                        fix_cage_geom=True,
                        sp_cage_subst=True,
                        heatplot=True,
                        arch_name='m2l4')


