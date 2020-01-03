from cgbind import Cage, Linker
from multiprocessing import Pool


# Define a function which will be called in parallel
def gen_cage(linker_name, linker_smiles):

    linker = Linker(name=linker_name, smiles=linker_smiles, arch_name='m2l4')
    cage = Cage(linker, metal_charge=2, metal='Pd')
    cage.print_xyzfile()

    return cage


# Define a dictionary of linker names and SMILES strings
linkers = {'L-1': 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
           'L-2': 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=N2)=CC=CN=C1'}

# Construct the cages in parallel with Python multiprocessing
with Pool(processes=2) as pool:
    results = [pool.apply_async(gen_cage, (name, smiles)) for name, smiles in list(linkers.items())]

    # Set no timeout, but it should be complete in 5 seconds!
    cages = [res.get(timeout=None) for res in results]


"""
To add substrates to the cavity see ex4
"""
