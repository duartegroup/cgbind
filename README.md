[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3550545.svg)](https://doi.org/10.5281/zenodo.3550545) [![Build Status](https://travis-ci.org/duartegroup/cgbind.svg?branch=master)](https://travis-ci.org/duartegroup/cgbind)
![alt text](cgbind/common/llogo.png)
### *C*a*g*e*Bind*ing: A tool to for automated metallocage binding affinity calculations
***

**cgbind** automates the process of generating and analysing metallocage structures
from just SMILES strings of the corresponding linkers. Binding affinity
calculations can be performed by input of both linker and substrate
SMILES. A graphical user interface (GUI) is available at
[cgbind.chem.ox.ac.uk](http://cgbind.chem.ox.ac.uk).

***

## Requirements
0. Python v. 3.7
1. All Python packages listed in requirements.txt 
2. [autode](https://github.com/duartegroup/autodE) latest
3. [ORCA](https://sites.google.com/site/orcainputlibrary/home) v. 4.1 (optional)
4. [XTB](https://github.com/grimme-lab/xtb) v. 6.2 (optional)

***

## Installation

For installation instructions see the [docs](https://duartegroup.github.io/cgbind/install.html).
If the requirements are already satisfied
```
python setup.py install
```

***

## Usage

**cgbind** constructs cages from [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)
strings, which can be generated directly from ChemDraw, as

![alt text](cgbind/common/SMILES_generation.png)


To generate a simple Pd<sub>2</sub>L<sub>4</sub> construct and print the geometry geometry
```python
from cgbind import Linker, Cage

linker = Linker(smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1', arch_name='m2l4')

cage = Cage(linker, metal='Pd')
cage.print_xyzfile()
```

See _examples/_ and the [docs](https://duartegroup.github.io/cgbind/examples.html)
for further examples.
