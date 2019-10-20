
class Config(object):
    #
    suppress_print = False
    #
    # Total number of cores available
    #
    n_cores = 1
    #
    # Supported codes: orca, xtb
    #
    code = None
    #
    # ------------------------------------------------------------------------------------------------
    # Available methods:
    #                     code = orca
    #                           0. pm3              : Semi-empirical PM3 method
    #                           1. xTB              : Tight-binding DFT method of Grimme
    #        DEFAULT            2. min_dft          : Minimal DFT method available in ORCA (PBE-D3BJ/def2-SVP,LooseOpt)
    #                           3. min_dft_normal   : Minimal DFT method available in ORCA (PBE-D3BJ/def2-SVP,NormalOpt)
    #                           4. dft              : PBE0-D3BJ/def2-SVP
    #
    #                     code = xtb
    #        DEFAULT            1. xtb      : Tight-binding method v 6.0 tested
    #
    #
    opt_method = None
    #
    #
    sp_keywords = ['SP', 'M062X', 'def2-TZVP', 'RIJCOSX', 'def2/J', 'SlowConv']
    sp_solvent = 'CH2Cl2'
    #
    # ------------------------------------------------------------------------------------------------
    # Paths:
    # The full path to an executable needs to be provided for any electronic structure method invoked
    #
    path_to_orca = None
    path_to_xtb = None
