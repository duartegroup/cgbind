import numpy as np
from cgbind.utils import work_in_tmp_dir
from cgbind.input_output import xyzfile_to_atoms
import re
from subprocess import Popen, DEVNULL, PIPE, STDOUT



class obabel():
    # TODO are keywords allowed
    # TODO Check if babel exists
    keywords = {"crit": 1e-6, 'ff': 'UFF', 'steps': 25000, 'rvdw': 10, 'rele': 10, 'freq': 10, 'log': True}
    energy = None

    def __init__(self, molecule=None):
        self.molecule = molecule

    def command_from_keywords(self, keywords):
        command = "obabel -imol2 temp.mol2 -oxyz -O out.xyz".split()
        for keyword in keywords:
            command.append(f'--{keyword:}')
            if not isinstance(keywords[keyword], bool):
                command.append(f'{keywords[keyword]}')
        return command

    def run_obabel(self, keywords):
        self.molecule.print_to_file(filename="temp.mol2")
        command = self.command_from_keywords(keywords)

        with open("output.txt", 'w') as output_file:
            process = Popen(command, stdout=DEVNULL, stderr=output_file)
            process.wait()

    @work_in_tmp_dir()
    def optimise(self, keywords={}):
        new_keywords = {**self.keywords, **keywords, **{'minimize': True}}
        self.run_obabel(new_keywords)
        atoms = xyzfile_to_atoms("out.xyz")
        self.molecule.set_atoms(atoms)

    @work_in_tmp_dir()
    def single_point(self, keywords={}):
        new_keywords = {**self.keywords, **keywords, **{'energy': True}}
        self.run_obabel(new_keywords)

        with open("output.txt", 'r') as File:
            text = File.read()
            if 'TOTAL ENERGY' in text:
                self.energy = float(re.search('TOTAL ENERGY = (\d+\.\d+) kJ', text).group(1))
            else:
                self.energy = None

    def get_energy(self):
        if self.energy is None:
            self.single_point()
        return self.energy

    def get_final_atoms(self):
        return self.molecule.atoms


ob = obabel()