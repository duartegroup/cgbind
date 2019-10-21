from cgbind.log import logger
from cgbind import add_substrate
from cgbind.input_output import print_output
from cgbind.geom import is_geom_reasonable
from cgbind import calculations
from cgbind.input_output import xyzs2xyzfile


class CageSubstrateComplex:

    def _reasonable_cage_substrate(self, cage, substrate):

        if cage is None or substrate is None:
            logger.error('Cannot build a cage-substrate complex either cage or substrate was None')
            return False

        attrs = [cage.charge, substrate.charge]
        if not all([attr is not None for attr in attrs]):
            logger.error('Cannot build a cage-substrate complex a required attribute was None')
            return False

        return True

    def print_xyzfile(self, force=False):
        if self.reasonable_geometry or force:
            xyzs2xyzfile(xyzs=self.xyzs, basename=self.name)

    def add_substrate(self):
        """
        Add a substrate to a cage.
        The binding mode will be determined by the number of heteroatoms etc. in the case of an M2L4 cage,
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

    def singlepoint(self, method, keywords, n_cores=1, max_core_mb=1000):
        return calculations.singlepoint(self, method, keywords, n_cores, max_core_mb)

    def optimise(self, method, keywords, n_cores=1, max_core_mb=1000, cartesian_constraints=None):
        return calculations.optimise(self, method, keywords, n_cores, max_core_mb, cartesian_constraints)

    def __init__(self, cage, substrate, solvent=None, mult=1):

        self.name = cage.name + '_' + substrate.name
        self.solvent = solvent
        self.mult = mult
        self.xyzs = None
        self.energy = None
        self.reasonable_geometry = False

        if not self._reasonable_cage_substrate(cage, substrate):
            return

        self.cage = cage
        self.substrate = substrate
        self.charge = cage.charge + substrate.charge

        self.add_substrate()
