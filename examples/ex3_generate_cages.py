# Import all the functions from cgbind
from cgbind import *

# Set the total available number of processing cores to 2 (default = 1)
Config.n_cores = 2

# Generate cages from L-1 and L-2 and don't optimise either of their geometries
gen_cages(linker_dict={
         'L-1'	: 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
         'L-2'	: 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=N2)=CC=CN=C1'},
          opt_linker=False,
          opt_cage=False)

"""
To add substrates to the cavity see ex4
"""
