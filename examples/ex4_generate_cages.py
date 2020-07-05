from cgbind import Cage, Linker, Config

# Use two cores. Both the RDKit conformer generation and components in the
# code are parallelised
Config.n_cores = 2


# Define a function to generate a cage
def gen_cage(linker_name, linker_smiles):

    linker = Linker(name=linker_name, smiles=linker_smiles, arch_name='m2l4')
    cage = Cage(linker, metal_charge=2, metal='Pd')
    cage.print_xyz_file()

    return cage


# Define a dictionary of linker names and SMILES strings
linkers = {'L-1': 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
           'L-2': 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=N2)=CC=CN=C1'}

# Construct the cages in serial
for name, smiles in linkers.items():
    gen_cage(linker_name=name, linker_smiles=smiles)

"""
To add substrates to the cavity see ex4
"""
