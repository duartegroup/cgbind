"""
Full list of functions available in cgbind, with all the function arguments set
"""

from cgbind import *

if __name__ == '__main__':

    Config.suppress_print = False
    Config.code = 'orca'
    Config.path_to_orca = '/Users/tom/opt/orca_4_1/orca'
    Config.opt_method = 'min_dft'
    Config.sp_keywords = ['SP', 'PBE', 'RI', 'def2-SVP']

    Config.path_to_mopac = '/Users/tom/opt/mopac/MOPAC2016.exe'
    Config.path_to_mopac_licence = '/Users/tom/opt/mopac/'

    Config.n_cores = 2

    Config.path_to_opt_struct = '/Users/tom/repos/cgbind/library/orca_opt_min_dft'

    if True:
        gen_cage('L-1', 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                 opt_linker=False,
                 opt_cage=False,
                 n_cores_pp=1)

    if True:
        gen_cages(linker_dict={
                 'L-1'	: 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                 'L-2'	: 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=N2)=CC=CN=C1',
                  },
                  opt_linker=True,
                  opt_cage=False
                  )

    if True:
        gen_cage_subst_complex('L-1', 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                               'q1', 'O=C1C=CC(C=C1)=O',
                               opt_cage=True,
                               opt_linker=False,
                               opt_substrate=False,
                               opt_cage_subst=False,
                               fix_cage_geom=True,
                               n_cores_pp=1)

    if True:
        gen_cage_subst_complexes(linker_dict={
                                'L-1'	: 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                                'L-2'	: 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=N2)=CC=CN=C1',
                                 },
                                 substrate_dict={'q1': 'O=C1C=CC(C=C1)=O'},
                                 opt_linker=False,
                                 opt_cage=False,
                                 opt_cage_subst=False,
                                 )

    if True:
        calc_binding_affinity('L-1', 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                              'q1', 'O=C1C=CC(C=C1)=O',
                              opt_cage=False,
                              opt_linker=False,
                              opt_substrate=True,
                              opt_cage_subst=False,
                              fix_cage_geom=False,
                              units_kcal_mol=True,
                              n_cores_pp=1)

    if True:
        calc_binding_affinities({'L-1'	: 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                                'L-2'	: 'C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=N2)=CC=CN=C1'},
                                {'quinone' : 'O=C1C=CC(C=C1)=O'},
                                substrate_charge=0, opt_linker=False, opt_cage=False,
                                opt_substrate=False, opt_cage_subst=False, linker_charge=0, metal_label='Pd',
                                metal_charge=2, fix_cage_geom=False, units_kcal_mol=False,
                                units_kj_mol=True, heatplot=True)
