# Import all the functions from cgbind
from cgbind import *

# Generate a cage from a linker name and SMILES string (which can be generated in Chemdraw using Edit/Copy As/SMILES)
gen_cage(linker_name='L-1',
         linker_smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1')

"""
To optimise either the linker before a cage is generated, or the full cage an electronic structure method needs to be
available and the path set in Config. See ex1
"""
