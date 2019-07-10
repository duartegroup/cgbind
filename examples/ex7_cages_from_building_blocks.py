from cgbind import *
from library.buildingblocks import *

Config.suppress_print = False
Config.code = 'mopac'
Config.path_to_mopac = '/Users/tom/opt/mopac/MOPAC2016.exe'
Config.path_to_mopac_licence = '/Users/tom/opt/mopac/'

Config.n_cores = 4

all_linkers = {}

# Loop over all possibilities for L components and build SMILES strings
for top_name, top_smiles in top_end_names_smiles.items():
    for link1_name, link1_smiles in top_link_names_smiles.items():
        for centre_name, centre_smiles in centre_names_smiles.items():
            for link2_name, link2_smiles in bottom_link_names_smiles.items():
                for end_name, end_smiles in bottom_end_names_smiles.items():
                    name = top_name + link1_name + centre_name + link2_name + end_name
                    smiles = make_linker_smiles(top_smiles, link1_smiles, centre_smiles, link2_smiles, end_smiles)

                    # Add all the unique combinations to the dictionary
                    if name not in all_linkers.keys():
                        all_linkers[name] = smiles

gen_cages(linker_dict=all_linkers,
          opt_linker=False,
          opt_cage=False,
          sp_cage=False
          )
