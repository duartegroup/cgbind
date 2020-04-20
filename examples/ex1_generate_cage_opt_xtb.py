from cgbind import Linker, Cage, XTB

# Generates the same linker as in ex0
linker = Linker(name='linker_opt', smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1', arch_name='m2l4')

# Optimise the linker with XTB using 2 cores
linker.optimise(method=XTB,
                keywords=None,
                n_cores=2)

# From the optimsied linker build the M2L4 cage with Pd(II) ions
cage = Cage(linker, metal='Pd', metal_charge='2')
cage.print_xyz_file()

"""
To optimise either the linker or the full cage at the DFT level see ex2
"""
