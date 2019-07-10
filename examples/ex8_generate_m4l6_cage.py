from cgbind import *

# Generate an M4L6 cage e.g. the cage presented in:
#               Caulder, D. L.; Powers, R. E.; Parac, T. N.; Raymond, K. N.. Angew. Chem., Int. Ed. 1998, 37, 1840−1843.

gen_cage(linker_name='Lx',
         linker_smiles='O=C(C1=C([O-])C([O-])=CC=C1)NC2=CC=CC3=C2C=CC=C3NC(C4=C([O-])C([O-])=CC=C4)=O',
         metal_label='Ga',
         metal_charge=3,
         arch='m4l6')
