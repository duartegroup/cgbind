from cgbind import input_output
import os


def test_io():

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

