from cgbind.linker import Linker
from cgbind.cage import Cage


def test_fragmentation_build():

    linker = Linker(name='linker',
                    # smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                    smiles='C1(/C=N/C2=CC=C(C=C2)C3=CC=C4C5=C(CC4=C3)C=C(C6=CC=C(C=C6)/N=C/C7=CC=CC=N7)C=C5)=CC=CC=N1',
                    # smiles='C1(C2=CC=C(C#CC3=CN=C(C4=NC=CC=C4)C=C3)C=N2)=NC=CC=C1',
                    arch_name='m4l6n',
                    use_fragment_conf=True)

    cage = Cage(linker, metal_charge=2, metal='Pd')
    cage.print_xyz_file()
