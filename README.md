[![pytest](https://github.com/duartegroup/cgbind/actions/workflows/pytest.yml/badge.svg)](https://github.com/duartegroup/cgbind/actions/workflows/pytest.yml) [![PyPI version](https://badge.fury.io/py/cgbind.svg)](https://badge.fury.io/py/cgbind)
![alt text](cgbind/common/llogo.png) 

**cgbind** automates the process of generating and analysing metallocage structures
from just SMILES strings of the corresponding linkers. Binding affinity
calculations can be performed with SMILES or 3D structures of linkers and the 
guest(s). A graphical user interface (GUI) is available at
[cgbind.chem.ox.ac.uk](http://cgbind.chem.ox.ac.uk).

## Installation

Install with:
```
conda install numpy autode rdkit scipy networkx Cython -c conda-forge
pip install cgbind
```
For detailed instructions see the [installation instructions](https://duartegroup.github.io/cgbind/install.html)

### Requirements
0. Python > v. 3.6
1. All Python packages listed in [requirements.txt](https://github.com/duartegroup/cgbind/blob/master/requirements.txt)
3. [ORCA](https://sites.google.com/site/orcainputlibrary/home) >v. 4.0 (optional)
4. [XTB](https://github.com/grimme-lab/xtb) >v. 6.2 (optional)

## Usage

**cgbind** constructs cages from [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)
strings, which can be generated directly from ChemDraw, as

![alt text](cgbind/common/SMILES_generation.png)


To generate a simple Pd<sub>2</sub>L<sub>4</sub> construct and print the geometry
```python
from cgbind import Linker, Cage

linker = Linker(smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1', arch_name='m2l4')

cage = Cage(linker, metal='Pd')
cage.print_to_file()
```

See [_examples/_](https://github.com/duartegroup/cgbind/tree/master/examples) and the 
[documentation](https://duartegroup.github.io/cgbind/examples.html) for further examples.


## Citation

If you find **cgbind** useful in your research please consider citing the paper:

T. A. Young, R. Gheorghe and F. Duarte, 
*cgbind: A Python Module and Web App for Automated Metallocage Construction and Host–Guest Characterization*.
J. Chem. Inf. Model. 2020, **60**, 7, 3546–3557. [DOI](https://doi.org/10.1021/acs.jcim.0c00519)