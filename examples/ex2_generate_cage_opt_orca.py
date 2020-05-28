from cgbind import Linker, Cage, orca

# Generates the same linker as in ex0
linker = Linker(name='linker_orca_opt', smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1', arch_name='m2l4')

# Optimise the linker with ORCA using 4 cores and fast DFT (PBE-D3BJ/def2-SVP)
linker.optimise(method=orca, n_cores=4)

# To specify a different set of keywords for a DFT optimisation at e.g.
# PBE0-D3BJ/def2-SVP use
#
# linker.optimise(method=orca,
#                 keywords=['LooseOpt', 'PBE0', 'D3BJ', 'def2-SVP'])

# From the optimsied linker build the M2L4 cage with Pd(II) ions
cage = Cage(linker, metal='Pd', metal_charge='2')
cage.print_xyz_file()

"""
To generate multiple cages in parallel see ex3
"""
