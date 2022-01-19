from cgbind import Linker, Cage

# Generate a simple cage
cage = Cage(linker=Linker(name='linker',
                          smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1',
                          arch_name='m2l4'),
            metal='Pd',
            metal_charge='2')

# Print the xyz file which can be loaded by e.g. Chimera and the ESP map over
# the VdW surface
cage.print_to_file()

# Generate and print a .cube file of the electrostatic potential using XTB
# partial atomic charges which provides qualitatively identical results to
# DFT and XTB integrals over the electron density
cage.print_esp_cube_file()
