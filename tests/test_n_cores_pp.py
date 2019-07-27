import os
from cgbind import *


def test_n_cores_pp():

    from cgbind.parallel import calc_n_cores_pp

    dict1, dict2 = {'test1': 'none', 'test2': 'none'}, {'test1': 'none', 'test2': 'none'}

    # Default total number of cores is 1
    assert calc_n_cores_pp(dict1, dict2) == 1

    Config.n_cores = 2
    assert calc_n_cores_pp(dict1) == 1
    assert calc_n_cores_pp(dict1, dict2) == 1

    Config.n_cores = 4
    assert calc_n_cores_pp(dict1) == 2
    assert calc_n_cores_pp(dict1, dict2) == 1
