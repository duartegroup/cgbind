from cgbind import Linker, Cage

# Initialise a cage from individual linkers to generate a M2L2L'2 cage

linker1 = Linker(name='L-1', smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1', arch_name='m2l4')
linker2 = Linker(name='L-2', smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=N2)=CC=CN=C1', arch_name='m2l4')

# Construct a metallocage from those linkers
cage = Cage(linkers=[linker1, linker2, linker1, linker2], metal='Pd', metal_charge=2)
cage.print_xyzfile()
