from cgbind.architectures import archs
import os


def test_arch():
    lib_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'cgbind', 'lib'))

    assert os.path.exists(lib_path)
    assert len(os.listdir(lib_path)) > 0
    assert 'm2l4.obj' in os.listdir(lib_path)

    assert any([arch.name == 'm2l4' for arch in archs])

    m2l4_arch = [arch for arch in archs if arch.name == 'm2l4'][0]
    assert m2l4_arch.n_metals == 2
    assert m2l4_arch.n_linkers == 4
