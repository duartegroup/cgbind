# Import the relevant objects
from cgbind import Linker, Cage

# Generate a linker name and SMILES string (which can be generated in Chemdraw
# using Edit/Copy As/SMILES)
linker = Linker(name='linker',
                smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                arch_name='m2l4')

# From the linker with a defined architecture build the cage with Pd(II) ions
cage = Cage(linker, metal='Pd', metal_charge='2')

# Print a .xyz file from the cage's xyzs
cage.print_to_file()


"""
To optimise a cage an electronic structure method needs to be available. 
See ex2 to optimise a cage using XTB

To see the progress of cgbind in your terminal run:

export CGBIND_LOG_LEVEL=INFO
"""
