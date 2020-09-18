from cgbind.linker import Linker
from cgbind.cage import Cage


def test_fragmentation():

    linker = Linker(name='linker',
                    smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                    arch_name='m2l4',
                    use_fragment_conf=True)

    cage = Cage(linker, metal_charge=2, metal='Pd')
    cage.print_xyz_file()
