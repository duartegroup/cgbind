from cgbind import Linker, Cage, CageSubstrateComplex, Substrate, xtb, Constants

# Calculate two binding affinities using XTB single points on generated
# structures. This will not be very accurate!

substrate = Substrate(name='quinone', smiles='O=C1C=CC(C=C1)=O')
substrate.singlepoint(method=xtb, n_cores=2)

# Make a list of linkers
linkers = [Linker(name='L-1', smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1', arch_name='m2l4'),
           Linker(name='L-2', smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=N2)=CC=CN=C1', arch_name='m2l4')]

# Compute binding affinities in serial
for linker in linkers:
    # Construct the cage
    cage = Cage(linker, metal='Pd', metal_charge=2)
    # Perform a single point energy evaluation at GFN-XTB
    cage.singlepoint(method=xtb, n_cores=2)

    # Construct the cage substrate complex
    cage_subt = CageSubstrateComplex(cage, substrate,
                                     energy_method='electrostatic_fast')
    # Perform a single point energy evaluation at GFN-XTB
    cage_subt.singlepoint(method=xtb, n_cores=2)

    binding_affinity = Constants.ha2kcalmol * (cage_subt.energy - cage.energy - substrate.energy)
    print('∆E_bind =', binding_affinity, 'kcal mol-1')
