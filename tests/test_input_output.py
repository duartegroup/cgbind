from cgbind import input_output
import os


def test_xyz_file():

    path_to_data = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))
    os.chdir(path_to_data)

    xyzs = input_output.xyzfile2xyzs(filename='cage_tmp.xyz')
    assert len(xyzs) == 146
    assert type(xyzs) == list
    assert type(xyzs[0]) == list
    assert type(xyzs[0][0]) == str
    assert type(xyzs[0][1]) == float

    input_output.xyzs2xyzfile(xyzs=xyzs, basename='tmp')
    assert os.path.exists('tmp.xyz')
    os.remove('tmp.xyz')

    # An xyz file should not be generated if..
    input_output.xyzs2xyzfile(xyzs=xyzs, basename=None)
    input_output.xyzs2xyzfile(xyzs=None, basename='tmp')
    input_output.xyzs2xyzfile(xyzs=xyzs, filename=None)
    assert not os.path.exists('tmp.xyz')


def test_mol_file():

    xyzs = input_output.molfile2xyzs(filename='methane.mol')
    assert len(xyzs) == 5
    assert type(xyzs[0]) == list
    assert type(xyzs[0][0]) == str
    assert type(xyzs[0][1]) == float


def test_mol2_file():

    xyzs = input_output.mol2file2xyzs(filename='methane.mol2')
    assert len(xyzs) == 5
    assert type(xyzs[0]) == list
    assert type(xyzs[0][0]) == str
    assert type(xyzs[0][1]) == float
