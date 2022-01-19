from cgbind import Linker, Cage

# Generate an M4L6 cage e.g. the cage presented in:
#  Caulder, D. L.; Powers, R. E.; Parac, T. N.; Raymond, K. N..
#  Angew. Chem., Int. Ed. 1998, 37, 1840−1843.

linker = Linker(name='linker_m4l6',
                smiles='O=C(C1=C([O-])C([O-])=CC=C1)NC2=CC=CC3=C2C=CC=C3NC(C4=C([O-])C([O-])=CC=C4)=O',
                arch_name='m4l6')

cage = Cage(linker, metal='Ga', metal_charge=3)
cage.print_to_file()
