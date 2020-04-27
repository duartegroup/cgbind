from cgbind import Cage, Linker, Substrate, CageSubstrateComplex

# Generate a cage.substrate complex
linker = Linker(name='L-1', smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1', arch_name='m2l4')
cage = Cage(linker, metal_charge=2, metal='Pd')

# Generate a substrate object
substrate = Substrate(name='benzoquione', smiles='O=C1C=CC(C=C1)=O')

# Build the cage.substrate complex from the cage and substrate objects
cs_complex = CageSubstrateComplex(cage=cage, substrate=substrate, energy_method='repulsion')
cs_complex.print_xyz_file()

"""
For a better binding mode use energy_method='electrostatic_fast'

To calculate binding affinities see ex5
"""
