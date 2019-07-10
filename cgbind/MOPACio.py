import os
from .log import logger
from .config import Config
from .constants import Constants
from subprocess import Popen


def singlepoint(xyzs, name, charge=0, mult=1):
    """
    Run a single point calculation using MOPAC
    :param xyzs: list of xyzs
    :param name: name of the species, will be appended with '_mopac_sp' to form the filename
    :param charge: total charge of the species
    :param mult: multiplicity of the species
    :return:
    """
    logger.info('Running MOPAC single point')

    energy = 0
    mop_filename = name + '_mopac_sp.mop'
    mop_out_filename = mop_filename.replace('.mop', '.out')

    gen_mopac_mop(mop_filename, xyzs, charge, mult, opt=False)
    if not Config.suppress_print:
        print("{:<30s}{:<50s}{:>10s}".format('Single point calculation of', name, 'Running'))
    mopac_out_lines = run_mopac(mop_filename, mop_out_filename)
    if not Config.suppress_print:
        print("{:<30s}{:<50s}{:>10s}".format('', name, 'Done'))

    for line in mopac_out_lines:
        if 'TOTAL ENERGY' in line:
            energy_ev = float(line.split()[-2])
            energy = energy_ev * Constants.eV2ha

    if energy == 0:
        logger.warning('Didn\'t find an energy in MOPAC output')

    return energy


def gen_mopac_mop(mop_filename, xyzs, charge, mult, opt=True, opt_atom_ids=None):
    logger.info('Generating a MOPAC input file')

    keywords = ['AUX', 'THREADS=1', 'CHARGE=' + str(charge)]
    if mult == 1:
        keywords.append('SINGLET')
    if mult == 2:
        keywords.append('DOUBLET')
    if mult == 3:
        keywords.append('TRIPLET')
    if Config.opt_method == 'pm6':
        keywords.append('PM6')
    elif Config.opt_method == 'pm7':
        keywords.append('PM7')
    else:
        keywords.append('PM7')
    if not opt:
        keywords.append('1SCF')

    with open(mop_filename, 'w') as mopac_mop:
        print(*keywords, '\nTitle \n', file=mopac_mop)
        if opt_atom_ids:
            for i in range(len(xyzs)):
                line = xyzs[i]
                if i in opt_atom_ids:
                    print(line[0], line[1], '1', line[2], '1', line[3], '1', sep='\t', file=mopac_mop)
                else:
                    print(line[0], line[1], '0', line[2], '0', line[3], '0', sep='\t', file=mopac_mop)
        else:
            [print(line[0], line[1], '1', line[2], '1', line[3], '1', sep='\t', file=mopac_mop) for line in xyzs]

    return 0


def run_mopac(mop_filename, out_filename):

    dev_null = open(os.devnull, 'w')
    mopac_run = Popen([Config.path_to_mopac, mop_filename],
                      stderr=dev_null,
                      env={'MOPAC_LICENSE': Config.path_to_mopac_licence})
    mopac_run.wait()

    return [line for line in open(out_filename, 'r')]


def get_mopac_xyzs_energy(out_lines):

    opt_converged, geom_section = False, False
    opt_xyzs = []
    energy = 0.0
    n_blank_lines = 0

    for line in out_lines:

        if 'SCF FIELD WAS ACHIEVED' in line:
            opt_converged = True

        if 'CARTESIAN COORDINATES' in line and opt_converged:
            geom_section = True

        if geom_section and len(line.split()) == 0:
            n_blank_lines += 1

        if geom_section and n_blank_lines == 2:
            geom_section = False

        if geom_section and len(line.split()) == 5:
            atom_label, x, y, z = line.split()[1:]
            opt_xyzs.append([atom_label, float(x), float(y), float(z)])

        if 'TOTAL ENERGY' in line:
            energy_ev = float(line.split()[-2])
            energy = Constants.eV2ha * energy_ev

    if len(opt_xyzs) == 0:
        logger.warning('Could not get optimised xyzs from MOPAC output')
    if energy == 0:
        logger.warning('Could not get optimised energy from MOPAC output')

    return opt_xyzs, energy
