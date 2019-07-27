from .log import logger
from . import m2l4
from .optimisation import opt_geom
from .single_point import singlepoint
from .input_output import print_output
from .geom import is_geom_reasonable
from .architectures import M2L4
from .architectures import M4L6


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
        self.xyzs = self.add_substrate_xyzs()

        if self.xyzs is not None:
            self.reasonable_geometry = is_geom_reasonable(self.xyzs)
        else:
            logger.warning('Cage-substrate xyzs are None')
            self.reasonable_geometry = False

        print_output('', self.substrate.name, 'Done')

    def add_substrate_xyzs(self):
        """
        Add the substrate to the cavity in a mode defined by the number of heteroatoms.

        :param substrate:
        :return:
        """
        logger.info('Adding substrate xyzs to cage')

        xyzs = None

        if self.cage.arch == M2L4:
            if self.substrate.x_x_atom_ids is not None:
                xyzs = m2l4.add_substrate_x_x(self.cage, self.substrate)
            elif self.substrate.x_atom_ids:
                xyzs = m2l4.add_substrate_x(self.cage, self.substrate)
            else:
                xyzs = m2l4.add_substrate_com(self.cage, self.substrate)

        if self.cage.arch == M4L6:
            # TODO add substrate with M4L6 cage
            logger.critical('Adding a substrate to an M4L6 cage. NOT IMPLEMENTED YET')
            exit()

        return xyzs

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