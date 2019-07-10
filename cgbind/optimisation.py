from .input_output import xyzs2xyzfile
from .ORCAio import *
from .XTBio import *
from .MOPACio import *


def opt_geom(xyzs, name='tmp', charge=0, mult=1, opt_atom_ids=None, n_cores=1):
    """
    Optimise the geometry of a generated species using an electronic structure code specified in config.py

    :param xyzs: List of xyzs
    :param name: Filename to use for the job
    :param charge: Total charge on the species
    :param mult: Multiplicity of the species
    :param opt_atom_ids: Atom numbers to optimise, keeping all others fixed. If None a full optimisation is performed
    :param n_cores: Number of cores to use for an ORCA optimisation
    :return: The new xyzs and the energy in Hartrees
    """

    energy = 0.0
    if xyzs is None:
        logger.error('Have no xyzs. Skipping optimisation')
        return None, energy

    if not Config.suppress_print:
        print("{:<30s}{:<50s}{:>10s}".format('Optimisation of', name, 'Running'))

    if Config.code == 'orca':
        if Config.path_to_orca is None:
            logger.error('path_to_orca needs to be set for an ORCA optimisation. Skipping the optimisation')
            if not Config.suppress_print:
                print("{:<30s}{:<50s}{:>10s}".format('', '', 'Failed'), '\n')
            return xyzs, energy

        inp_filename = name + '_orca_opt.inp'
        out_filename = inp_filename.replace('.inp', '.out')

        gen_orca_inp(inp_filename, xyzs, charge, mult, opt=True, opt_atom_ids=opt_atom_ids, n_cores=n_cores)
        orca_output_file_lines = run_orca(inp_filename, out_filename)
        opt_xyzs, energy = get_orca_xyzs_energy(out_lines=orca_output_file_lines)

    elif Config.code == 'xtb':
        if Config.path_to_xtb is None:
            logger.error('path_to_xtb needs to be set for a XTB optimisation. Skipping the optimisation')
            if not Config.suppress_print:
                print("{:<30s}{:<50s}{:>10s}".format('', '', 'Failed'), '\n')
            return xyzs, energy

        xtb_filename = name + '_xtb_opt'
        xtb_input_file = xyzs2xyzfile(xyzs, xtb_filename + '.xyz')
        xtb_output_lines = run_xtb(xtb_input_file, opt=True, charge=charge)
        opt_xyzs, energy = get_XTB_xyzs_energy(out_lines=xtb_output_lines)

    elif Config.code == 'mopac':
        if Config.path_to_mopac is None or Config.path_to_mopac_licence is None:
            logger.error('path_to_mopac and path_to_mopac_licence need to be set. Skipping the optimisation')
            if not Config.suppress_print:
                print("{:<30s}{:<50s}{:>10s}".format('', '', 'Failed'), '\n')
            return xyzs, energy

        mop_filename = name + '_mopac_opt.mop'
        mop_out_filename = mop_filename.replace('.mop', '.out')

        gen_mopac_mop(mop_filename, xyzs, charge, mult, opt=True, opt_atom_ids=opt_atom_ids)
        mopac_output_lines = run_mopac(mop_filename, mop_out_filename)
        opt_xyzs, energy = get_mopac_xyzs_energy(out_lines=mopac_output_lines)

    elif Config.code is None:
        logger.error('An optimisation was called but no code was set. Skipping the optimisation')
        if not Config.suppress_print:
            print("{:<30s}{:<50s}{:>10s}".format('', '', 'Failed'), '\n')
        return xyzs, energy

    else:
        logger.error('No available code with name "{}". Skipping the optimisation'.format(Config.code))
        if not Config.suppress_print:
            print("{:<30s}{:<50s}{:>10s}".format('', '', 'Failed'), '\n')
        return xyzs, energy

    if len(opt_xyzs) == 0:
        logger.error('No xyz lines found. Returning the non-optimised xyz')
        return xyzs, energy

    if not Config.suppress_print:
        print("{:<30s}{:<50s}{:>10s}".format('', name, 'Done'), '\n')

    return opt_xyzs, energy
