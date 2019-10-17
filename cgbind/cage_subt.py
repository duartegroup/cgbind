from cgbind.log import logger
from cgbind import add_substrate
from cgbind.optimisation import opt_geom
from cgbind.single_point import singlepoint
from cgbind.input_output import print_output
from cgbind.geom import is_geom_reasonable


class CageSubstrateComplex(object):

    def add_substrate(self):
        """
        Add a substrate to a cage.
        The binding mode will be determined by the number of heteroatoms etc. in the case of an M2L4 cage,
        :param substrate: Substrate object
        :return:
        """

        if self.cage.xyzs is None or self.substrate.xyzs is None:
            logger.error("Cage and/or substrate has no xyzs. Can't add a substrate")
            return

        print_output('Addition of', self.substrate.name, 'Running')
        logger.info('Adding the substrate to the center of the cage defined by the COM')
        self.xyzs = add_substrate.add_substrate_com(self.cage, self.substrate)

        if self.xyzs is not None:
            self.reasonable_geometry = is_geom_reasonable(self.xyzs)
        else:
            logger.error('Cage-substrate xyzs are None')
            self.reasonable_geometry = False

        print_output('', self.substrate.name, 'Done')

    def optimise(self, opt_atom_ids=None, n_cores=1):
        """
        Like optimise, but with a cage.substrate species
        :param opt_atom_ids: Atom ids to optimise, possible to just optimise the substrate within the cage,
        providing a considerable computational speedup
        :param n_cores: Number of cores to use for the optimisation
        :return:
        """
        logger.info('Optimising a cage-substrate complex')
        self.xyzs, self.energy = opt_geom(self.xyzs, self.name, charge=self.charge, opt_atom_ids=opt_atom_ids,
                                          n_cores=n_cores)

    def singlepoint(self, n_cores=1):
        self.energy = singlepoint(self, n_cores)

    def __init__(self, cage, substrate):

        self.cage = cage
        self.substrate = substrate
        self.xyzs = None
        self.reasonable_geometry = False
        self.name = cage.name + '_' + substrate.name
        self.charge = cage.charge + substrate.charge

        self.energy = None
        self.add_substrate()
