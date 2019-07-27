# Installation
## Windows
***
On Windows **cgbind** can be installed with the following steps: 

1.First download [git](https://git-scm.com/download/win) and follow the instructions with all the defaults

2.Download the Windows Minconda installer from [here](https://docs.conda.io/en/latest/miniconda.html)

3.Follow the instructions to install Miniconda. Make sure to tick the **Add Anaconda to my PATH environment variable** 
box

4.Launch _Git bash_ from the Start Menu 

5.Change directories to your Downloads folder by typing
```
cd Downloads/
```
and hitting enter.

6.Copy this directory wih
```
git clone https://github.com/duartegroup/cgbind.git
```

7.Change directories into _cgbind/_
```
cd cgbind/
```

8.Run 
```
conda config --append channels conda-forge && conda env create --file requirements.txt --name cgbind-env
```
which will take a few minute to download then install all the Python dependencies to a new virtual environment 
(_cgbind-env_)

9.Run
```
conda init bash
```
and reopen _Git bash_

10.Activate the virtual Python environment with
```
conda activate cgbind-env
```

7.Change directories back to _cgbind/_
```
cd Downloads/cgbind/
```

11.Finally, **cgbind** installed with
```
python setup.py install
```


## Mac
***
The easiest way to install *cgbind* on Mac is to use [Miniconda](https://docs.conda.io/en/latest/miniconda.html) 
to satisfy the dependencies automatically. See [here](https://vimeo.com/347275041) for an example fresh installation 
walk-though. A step-by step install:

1.Open up a terminal window and type 
```
cd Downloads/
```
and hit enter to change directories to your downloads folder

2.Run
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```
to grab the latest Miniconda version

3.Run
```
bash Miniconda3-latest-MacOSX-x86_64.sh
```
and follow the instructions (keep hitting enter and typing yes where appropriate) to install Miniconda. 

4.To install **cgbind** we'll need the files in this repository. Download them with git 
```
git clone https://github.com/duartegroup/cgbind.git
```
5.Open a **new** terminal window change directories to _cgbind/_
```
cd Downloads/cgbind/
```

6.Install the dependencies with
```
conda config --append channels conda-forge && conda env create --file requirements.txt --name cgbind-env
```
7.Activate the virtual Python environment with
```
conda activate cgbind-env
```
8.Finally, **cgbind** installed with
```
python setup.py install
```

Everything should now work! 

