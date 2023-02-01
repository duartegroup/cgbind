from typing import List, Optional
import numpy as np

from cgbind.log import logger
from cgbind.config import Config
from cgbind.exceptions import RequiresAutodE


def _cgbind_mol_to_autode(
        cgbind_mol: "cgbind.molecule.BaseStruct"
) -> "autode.Molecule":
    """
    Convert a cgbind Moleucle/Cage/BaseStruct object into an
    autodE Species that can be dispatched directly to autodE

    :param cgbind_mol: The cgbind BaseStruct object
    :return: (autode.Molecule)
    """
    # todo should it use autode.Species or autode.Molecule?
    from cgbind.molecule import BaseStruct
    try:
        from autode.atoms import Atom
        from autode import Molecule
    except ModuleNotFoundError:
        logger.error('autode not found. Calculations are not available')
        raise RequiresAutodE  # should other errors be removed?

    assert isinstance(cgbind_mol, BaseStruct)

    atom_list = []
    for atom in cgbind_mol.atoms:
        coords = list(atom.coord)
        assert len(coords) == 3
        ade_atom = Atom(atom.label, *coords)
        atom_list.append(ade_atom)

    ade_mol = Molecule(name=cgbind_mol.name, atoms=atom_list,
                       charge=cgbind_mol.charge, mult=cgbind_mol.mult,
                       solvent_name=cgbind_mol.solvent)

    if hasattr(cgbind_mol, "conformers"):
        if cgbind_mol.conformers is not None:
            for conf in cgbind_mol.conformers:
                ade_conf = _conformer_to_autode(conf, ade_mol)
                ade_mol.conformers.append(ade_conf)

    return ade_mol


def _conformer_to_autode(cgbind_conf: "cgbind.molecule.BaseStruct",
                         parent_ade_mol: "autode.Species"):
    """
    Take a conformer from a cgbind Molecule and convert it to
    autodE

    :param cgbind_conf: (cgbind.molecule.BaseStruct) a conformer
    :param parent_ade_mol: (autode.Species) the parent species
    :return:
    """
    # conformer cannot have more conformers!
    assert not hasattr(cgbind_conf, "conformers")
    try:
        from autode.conformers.conformer import Conformer
    except ModuleNotFoundError:
        raise RequiresAutodE

    ade_conf = Conformer(species=parent_ade_mol, name=cgbind_conf.name)
    converted_conf = _cgbind_mol_to_autode(cgbind_conf)
    ade_conf.atoms = converted_conf.atoms

    return ade_conf


def optimise(molecule, method, keywords, n_cores=None, cartesian_constraints=None):
    """
    Optimise a molecule

    :param molecule: (cgbind.molecule.BaseStruct)
    :param method: (autode.ElectronicStructureMethod)
    :param keywords: (list(str)) Keywords to use for the electronic structure c
                                 alculation e.g. ['Opt', 'PBE', 'def2-SVP']
    :param n_cores: (int) Number of cores to use
    :param cartesian_constraints: (list(int)) List of atom ids to constrain
    :return:
    """
    logger.info('Running an optimisation calculation')

    n_cores = Config.n_cores if n_cores is None else int(n_cores)
    ade_mol = _cgbind_mol_to_autode(molecule)

    try:
        from autode.calculation import Calculation
        from autode.wrappers.XTB import xtb
        from autode.wrappers.ORCA import orca
        from autode.wrappers.keywords import OptKeywords

    except ModuleNotFoundError:
        logger.error('autode not found. Calculations not available')
        raise RequiresAutodE

    if keywords is None:
        if method == orca:

            keywords = OptKeywords(['LooseOpt', 'PBE', 'D3BJ', 'def2-SVP'])
            logger.warning(f'No keywords were set for the optimisation but an'
                           f' ORCA calculation was requested. '
                           f'Using {str(keywords)}')

        elif method == xtb:
            keywords = xtb.keywords.opt

        else:
            logger.critical('No keywords were set for the optimisation '
                            'calculation')
            raise Exception

    else:
        # If the keywords are specified as a list convert them to a set of
        # OptKeywords, required for autodE
        if type(keywords) is list:
            keywords = OptKeywords(keywords)

    opt = Calculation(name=ade_mol.name + '_opt',
                      molecule=ade_mol,
                      method=method,
                      keywords=keywords,
                      n_cores=n_cores,
                      cartesian_constraints=cartesian_constraints)
    opt.run()
    molecule.energy = float(ade_mol.energy.to('Ha'))
    # todo is it sufficient to only update coordinates?
    opt_coords = np.array(ade_mol.coordinates.to('ang')).copy()
    assert opt_coords.shape == (molecule.n_atoms, 3)
    molecule.set_atoms(coords=opt_coords)

    return None


def singlepoint(molecule, method, keywords, n_cores=None):
    """
    Run a single point energy evaluation on a molecule

    :param molecule: (object)
    :param method: (autode.ElectronicStructureMethod)
    :param keywords: (list(str)) Keywords to use for the electronic structure
    calculation e.g. ['Opt', 'PBE', 'def2-SVP']
    :param n_cores: (int) Number of cores to use
    :return:
    """
    logger.info('Running single point calculation')

    n_cores = Config.n_cores if n_cores is None else int(n_cores)
    ade_mol = _cgbind_mol_to_autode(molecule)

    try:
        from autode.calculation import Calculation
        from autode.wrappers.XTB import xtb
        from autode.wrappers.ORCA import orca
        from autode.wrappers.keywords import SinglePointKeywords

    except ModuleNotFoundError:
        logger.error('autode not found. Calculations not available')
        raise RequiresAutodE

    if keywords is None:
        if method == orca:
            keywords = SinglePointKeywords(['SP', 'M062X', 'def2-TZVP',
                                            'RIJCOSX', 'def2/J', 'SlowConv'])

            logger.warning('No keywords were set for the single point but an '
                           'ORCA calculation was requested. '
                           f'Using {str(keywords)}')

        elif method == xtb:
            keywords = xtb.keywords.sp

        else:
            logger.critical('No keywords were set for the single-point '
                            'calculation')
            raise Exception

    else:
        # If the keywords are specified as a list convert them to a set of
        # OptKeywords, required for autodE
        if type(keywords) is list:
            keywords = SinglePointKeywords(keywords)

    sp = Calculation(name=ade_mol.name + '_sp', molecule=ade_mol,
                     method=method,
                     keywords=keywords,
                     n_cores=n_cores)
    sp.run()
    molecule.energy = float(ade_mol.energy.to('Ha'))

    return None


def get_charges(molecule) -> Optional[List[float]]:
    """
    Get the partial atomic charges with XTB (tested with v. 6.2) will generate
    then trash a temporary directory

    :return: (list(float)) a list of partial charges
    """
    logger.info('Getting charges')
    ade_mol = _cgbind_mol_to_autode(molecule)

    try:
        from autode.calculations import Calculation
        from autode.wrappers.XTB import xtb
        from autode.exceptions import MethodUnavailable
        from autode.wrappers.keywords import SinglePointKeywords

    except ModuleNotFoundError:
        logger.error('autode not found. Calculations not available')
        raise RequiresAutodE

    # Run the calculation
    try:
        xtb_sp = Calculation(name=ade_mol.name + '_xtb_sp', molecule=ade_mol,
                             method=xtb,
                             n_cores=1,
                             keywords=xtb.keywords.sp)
        xtb_sp.run()

        charges = [float(chrg) for chrg in ade_mol.partial_charges()]

    except MethodUnavailable:
        logger.error('Could not calculate without an XTB install')
        return None

    if len(charges) == molecule.n_atoms:
        return charges

    else:
        logger.error('XTB failed to generate charges')
        return None
