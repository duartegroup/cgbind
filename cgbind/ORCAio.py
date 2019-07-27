import os
from .log import logger
from .input_output import print_output
from .config import Config
from subprocess import Popen


def get_single_point_energy(xyzs, name, keywords, solvent=None, charge=0, mult=1, n_cores=1):
    logger.info('Running an ORCA single point of {}'.format(name))

    if not Config.code == 'orca':
        logger.warning('Code was set at "{}", but only ORCA single points are supported, '
                       'Using ORCA anyway'.format(Config.code))

    energy = 0
    inp_filename = name + '_orca_sp.inp'
    out_filename = inp_filename.replace('.inp', '.out')

    if ('PAL' + str(n_cores)) not in keywords and n_cores > 1:
        if n_cores > 8:
            logger.warning("Number of ORCA cores requested was {}. Reducing to 8".format(n_cores))
            keywords.append('PAL8')
        else:
            keywords.append('PAL' + str(n_cores))

    if solvent and 'CPCM' not in keywords:
        keywords.append('CPCM')

    xys2inp(xyzs, keywords, inp_filename, charge, mult, solvent)
    if Config.path_to_orca is None:
        logger.error('path_to_orca needs to be set for an ORCA single point. Skipping')
        return

    print_output('Single point calculation of', name, 'Running')
    orca_out_lines = run_orca(inp_filename, out_filename)
    print_output(' ', name, 'Done')
    for line in orca_out_lines:
        if 'FINAL SINGLE POINT ENERGY' in line:
            energy = float(line.split()[4])

    return energy


def xys2inp(xyzs, keywords, filename, charge=0, mult=1, solvent=None, opt_atom_ids=None):

    if filename.endswith('.inp'):
        with open(filename, 'w') as inp_file:
            print('!', *keywords, file=inp_file)
            if solvent:
                print('%cpcm\n smd true\n SMDsolvent \"' + solvent + '\"\n end', file=inp_file)
            print('%scf \nmaxiter 250 \nend', file=inp_file)
            print('% maxcore 4000', file=inp_file)
            if opt_atom_ids:
                print('%geom Constraints', file=inp_file)
                for i in opt_atom_ids:
                    print('{ C', i, 'C }', file=inp_file)
                print('end\n', 'invertConstraints true\n', 'end', sep='', file=inp_file)
            print('*xyz', charge, mult, file=inp_file)
            [print(*line, sep='\t', file=inp_file) for line in xyzs]
            print('*', file=inp_file)

    return 0


def gen_orca_inp(inp_filename, xyzs, charge, mult, opt=True, opt_atom_ids=None, n_cores=1):
    logger.info('Generating {}'.format(inp_filename))

    if Config.opt_method == 'xtb':
        keywords = ['LooseOpt', 'XTB']
    elif Config.opt_method == 'min_dft':
        keywords = ['LooseOpt', 'PBE', 'D3BJ', 'RI', 'def2-SVP', 'SlowConv']
    elif Config.opt_method == 'min_dft_normal':
        keywords = ['Opt', 'PBE', 'RI', 'D3BJ', 'def2-SVP', 'SlowConv']
    elif Config.opt_method == 'dft':
        keywords = ['Opt', 'PBE0', 'RIJCOSX', 'D3BJ', 'def2-SVP']
    else:
        logger.warning('No opt method set. Defaulting to PBE-D3BJ/def2-SVP/LooseOpt')
        keywords = ['LooseOpt', 'PBE', 'D3BJ', 'RI', 'def2-SVP', 'SlowConv']

    if not opt:
        if 'LooseOpt' in keywords:
            keywords.remove('LooseOpt')
        if 'Opt' in keywords:
            keywords.remove('Opt')

    if 1 < n_cores <= 8:
        keywords.append('PAL' + str(n_cores))
    if n_cores > 8:                                 # ORCA calculations are capped at using 8 cores
        logger.warning("Number of ORCA cores requested was {}. Reducing to 8".format(n_cores))
        keywords.append('PAL8')

    xys2inp(xyzs, keywords, inp_filename, charge, mult, solvent=None, opt_atom_ids=opt_atom_ids)

    return 0


def run_orca(inp_filename, out_filename):
    logger.info('Running ORCA calculation from {}'.format(inp_filename))

    orca_done, orca_terminated_normally = False, False

    if os.path.exists(out_filename):
        logger.info('Out file already exists. Checking to see whether it completed')
        orca_done = True

    if orca_done:
        out_lines = [line for line in open(out_filename, 'r', encoding="utf-8")]
        for line in reversed(out_lines):
            if 'ORCA TERMINATED NORMALLY' in line:
                orca_terminated_normally = True
                if not Config.suppress_print:
                    print("{:<30s}{:<50s}{:>10s}".format('Found ORCA run done for', inp_filename, 'Skipping'))
                break

    if not orca_terminated_normally:
        with open(out_filename, 'w') as orca_out:
            orca_run = Popen([Config.path_to_orca, inp_filename], stdout=orca_out)
        orca_run.wait()

    return [line for line in open(out_filename, 'r', encoding="utf-8")]


def get_orca_xyzs_energy(out_lines):
    logger.info('Getting xyzs and final energy from an ORCA output file')

    opt_converged, geom_section = False, False
    opt_xyzs = []
    energy = 0.0

    for line in out_lines:

        if 'THE OPTIMIZATION HAS CONVERGED' in line:
            opt_converged = True
        if 'CARTESIAN COORDINATES' in line and opt_converged:
            geom_section = True

        if geom_section and len(line.split()) == 0:
            geom_section = False

        if geom_section and len(line.split()) == 4:
            atom_label, x, y, z = line.split()
            opt_xyzs.append([atom_label, float(x), float(y), float(z)])

        if 'FINAL SINGLE POINT ENERGY' in line:
            energy = float(line.split()[4])             # e.g. line = 'FINAL SINGLE POINT ENERGY     -4143.815610365798'

    if len(opt_xyzs) == 0:
        logger.warning('Could not find any xyzs in ORCA output')
    if energy == 0:
        logger.warning('Could not find energy in ORCA output')

    return opt_xyzs, energy
