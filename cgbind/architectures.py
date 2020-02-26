import os
import pickle
from cgbind.log import logger


class Arch:

    def __init__(self, name, n_metals, n_linkers):
        """
        Architecture class holding the name (self.name), self.n_metals and self.n_linkers

        :ivar self.name: (str)
        :ivar self.n_metals: (int)
        :ivar self.n_linkers: (int)

        :param name: (str) Name of the architecture
        :param n_metals: (int) Number of metals in the architecture. e.g. 2 for M2L4
        :param n_linkers: (int) Number of linkers in the architecture. e.g. 4 for M2L4
        """
        self.name = name
        self.n_metals = n_metals
        self.n_linkers = n_linkers


archs = []

# Populate the available architectures from lib/
folder_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib')

logger.info('Populating available architectures')
for filename in os.listdir(folder_path):
    if filename.endswith('.obj'):
        pickled_file = open(os.path.join(folder_path, filename), 'rb')
        template = pickle.load(pickled_file)

        archs.append(Arch(name=template.arch_name,
                          n_metals=template.n_metals,
                          n_linkers=template.n_linkers))

logger.info(f'Have {len(archs)} architectures: {[arch.name for arch in archs]}')
