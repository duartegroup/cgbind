[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3930811.svg)](https://doi.org/10.5281/zenodo.3930811) [![Build Status](https://travis-ci.org/duartegroup/cgbind.svg?branch=master)](https://travis-ci.org/duartegroup/cgbind) [![PyPI version](https://badge.fury.io/py/cgbind.svg)](https://badge.fury.io/py/cgbind)
![alt text](cgbind/common/llogo.png) 

**cgbind** automates the process of generating and analysing metallocage structures
from just SMILES strings of the corresponding linkers. Binding affinity
calculations can be performed with SMILES or 3D structures of linkers and the 
guest(s). A graphical user interface (GUI) is available at
[cgbind.chem.ox.ac.uk](http://cgbind.chem.ox.ac.uk).

## Installation

For detailed instructions see the [installation instructions](https://duartegroup.github.io/cgbind/install.html).
If the requirements are already satisfied
```
pip install cgbind
```

### Requirements
0. Python v. 3.7
1. All Python packages listed in requirements.txt 
2. [autode](https://github.com/duartegroup/autodE) latest (optional)
3. [ORCA](https://sites.google.com/site/orcainputlibrary/home) v. 4.1 (optional)
4. [XTB](https://github.com/grimme-lab/xtb) v. 6.2 (optional)

The Python dependencies are most easily satisfied using a conda
([anaconda](https://www.anaconda.com/distribution)/[miniconda](https://docs.conda.io/en/latest/miniconda.html))
installation by running

```
conda config --append channels conda-forge
conda install --file requirements.txt
```


## Usage

**cgbind** constructs cages from [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)
strings, which can be generated directly from ChemDraw, as

![alt text](cgbind/common/SMILES_generation.png)


To generate a simple Pd<sub>2</sub>L<sub>4</sub> construct and print the geometry
```python
from cgbind import Linker, Cage

linker = Linker(smiles='C1(C#CC2=CC=CC(C#CC3=CC=CN=C3)=C2)=CC=CN=C1', arch_name='m2l4')

cage = Cage(linker, metal='Pd')
cage.print_xyz_file()
```

See _examples/_ and the [docs](https://duartegroup.github.io/cgbind/examples.html)
for further examples.


## Citation

If you find **cgbind** useful in your research please consider citing the paper:

T. A. Young, R. Gheorghe and F. Duarte, 
*cgbind: A Python Module and Web App for Automated Metallocage Construction and Host–Guest Characterization*.
J. Chem. Inf. Model. 2020, **60**, 7, 3546–3557. [DOI](https://doi.org/10.1021/acs.jcim.0c00519)