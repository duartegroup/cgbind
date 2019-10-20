from cgbind import *

Config.code = 'xtb'

Config.n_cores = 4

calc_binding_affinities({'L-1': 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                         'L-2': 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=N2)=CC=CN=C1'},
                        {'quinone': 'O=C1C=CC(C=C1)=O',
                         'dimethylquinone': 'O=C(C=C1C)C(C)=CC1=O'},
                        opt_linker=True,
                        opt_cage=True,
                        opt_substrate=True,
                        opt_cage_subst=True,
                        metal_label='Pd',
                        metal_charge=2,
                        fix_cage_geom=False,
                        sp_cage=False,
                        sp_substrate=False,
                        sp_cage_subst=False,
                        units_kcal_mol=False,
                        units_kj_mol=True,
                        heatplot=True)
