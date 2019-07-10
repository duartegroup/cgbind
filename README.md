![alt text](cgbind/common/llogo.png)

### *C*a*g*e*Bind*ing: A tool to for automated binding affinity calculations

Constructing metallocage structures 


## Requirements
0. Python v. 3.7
1. [rdkit](https://github.com/rdkit/rdkit)
2. matplotlib
3. ORCA v. 4.1 (optional)
4. xTB (optional)
5. MOPAC (optional)

## Installation

If the requirements are already satisfied
```
python setup.py install
```

Alternately for a one line fresh install with miniconda on Mac OSX
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh && bash Miniconda3-latest-MacOSX-x86_64.sh && source ~/.bash_profile && conda install -c rdkit rdkit && git clone https://tyoung31@bitbucket.org/duartegroup/cgbind.git && cd cgbind && python setup.py install
```


## Usage
To load as a module:
```
import cgbind
```

See _examples/_ for how to use _cgbind_
