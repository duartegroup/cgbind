from cgbind import Linker, Cage
from multiprocessing import Pool


# Function to generate a cage and print the associated xyz file
def gen_cage(linker_name, linker_smiles):

    linker = Linker(name=linker_name, smiles=linker_smiles, arch_name='m2l4')
    cage = Cage(linker, metal_charge=2, metal='Pd')
    cage.print_xyz_file()

    return cage


# Dictionaries of the building blocks where the M2L4 linker looks like:
#
#  top_end
#  |
#  top link
#  |
#  centre
#  |
#  bottom link
#  |
#  bottom end

top_end_names_smiles = {'3-pyridine': 'C1=CC=NC=C1%99'}

top_link_names_smiles = {'alkyne': '.C%99#C%98.'}

centre_names_smiles = {'2-6-pyridine'         : 'C1%98=NC%97=CC=C1',
                       'carbonate'            : 'O=C(O%98)O%97',
                       '1-3-phenyl'           : 'C1%98=CC=CC%97=C1',
                       '1-3-phenyl-2-methyl'  : 'C1%98=CC=CC%97=C1C',
                       'methylene'            : 'C%98%97'}

bottom_link_names_smiles = {'alkyne': '.C%97#C%96.'}

bottom_end_names_smiles = {'3-pyridine': 'C1%96=CC=CN=C1'}

# Empty dictionary to store the linker names and SMILES strings
all_linkers = {}

# Loop over all possibilities for L components and build SMILES strings
for top_name, top_smiles in top_end_names_smiles.items():
    for link1_name, link1_smiles in top_link_names_smiles.items():
        for centre_name, centre_smiles in centre_names_smiles.items():
            for link2_name, link2_smiles in bottom_link_names_smiles.items():
                for end_name, end_smiles in bottom_end_names_smiles.items():
                    name = top_name + link1_name + centre_name + link2_name + end_name
                    smiles = top_smiles + link1_smiles + centre_smiles + link2_smiles + end_smiles

                    # Add all the unique combinations to the dictionary
                    if name not in all_linkers.keys():
                        all_linkers[name] = smiles


# Parallelise over 4 cores
with Pool(processes=4) as pool:
    # Apply asynchronously because the cages are independent
    results = [pool.apply_async(gen_cage, (name, smiles)) for name, smiles in list(all_linkers.items())]

    cages = [res.get(timeout=None) for res in results]
