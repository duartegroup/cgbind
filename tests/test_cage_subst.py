from cgbind import cage
from cgbind.linker import Linker
from cgbind.substrate import Substrate
from cgbind import cage_subt
import os

here = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()


def test_cage_subst_gen():
    os.chdir(os.path.join(here, 'data'))

    linker = Linker(arch_name='m2l4', filename='ch_linker.xyz', name='tmp')
    c = cage.Cage(linker=linker, metal='Pd', metal_charge=2)
    s = Substrate(smiles='C', name='methane', n_confs=1)

    cs = cage_subt.CageSubstrateComplex(cage=c, substrate=s)

    assert cs.reasonable_geometry is True
    assert cs.charge == 4

    cs.print_xyzfile()
    assert os.path.exists('cage_tmp_methane.xyz')
    os.remove('cage_tmp_methane.xyz')

    os.chdir(cwd)
