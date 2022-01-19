"""
Due to either the conformational freedom of the linker or the imperfect
SMILES -> 3D geometry it may sometimes be more reliable to initialise cages
from 3D structures. For example, to generate the M4L6 cage used in
Angew. Chem. Int. Ed. 2008, 47, 8297 from a .xyz of the constituent linker
"""
from cgbind import Linker, Cage

linker = Linker(name='m4l6_linker',
                filename='m4l6_linker.xyz',
                arch_name='m4l6n')

cage = Cage(linker, metal='Fe', metal_charge='2')
cage.print_to_file()
