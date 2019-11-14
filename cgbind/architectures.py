import os
import pickle
from cgbind.log import logger


class Arch:

    def __init__(self, name, n_metals, n_linkers):

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

logger.info(f'Have {len(archs)} architectures')
