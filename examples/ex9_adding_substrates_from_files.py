from cgbind import Cage, Substrate, Linker, CageSubstrateComplex

# Generate an M4L6 bipyridyl Fe(III) based metallocage
linker = Linker(name='linker',
                smiles='C1(C2=CC=C(C#CC3=CN=C(C4=NC=CC=C4)C=C3)C=N2)=NC=CC=C1',
                arch_name='m4l6n')

cage = Cage(linker, metal='Fe', metal_charge='3')

# Initialise a substrate from the .xyz file, setting the charge and spin
# multiplicity
substrate = Substrate(name='PhO-', filename='phenoxide.xyz', charge=-1, mult=1)

cs = CageSubstrateComplex(cage, substrate)
cs.print_to_file()
