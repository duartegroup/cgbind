import os
from cgbind.log import logger
from cgbind.config import Config
from cgbind.constants import Constants
from subprocess import Popen


def run_xtb(xyz_filename, opt=True, charge=0, n_cores=1):
    logger.info('Running XTB')

    xtb_out_filename = xyz_filename.replace('.xyz', '.out')

    if os.path.exists(xtb_out_filename):
        out_lines = [line for line in open(xtb_out_filename, 'r')]
        for line in out_lines:
            if ' * finished run' in line:
                if not Config.suppress_print:
                    print("{:<30s}{:<50s}{:>10s}".format('Found XTB run done for', xyz_filename, 'Skipping'))
                return out_lines

    with open(xtb_out_filename, 'w') as xtb_out:
        if opt:
            params = [Config.path_to_xtb, xyz_filename, '--opt', '-c', str(charge)]
        else:
            params = [Config.path_to_xtb, xyz_filename, '-c ', str(charge)]

        os.environ['OMP_NUM_THREADS'] = str(n_cores)
        xtb_run = Popen(params, stdout=xtb_out, stderr=open(os.devnull, 'w'))
    xtb_run.wait()

    return [line for line in open(xtb_out_filename, 'r')]


def get_XTB_xyzs_energy(out_lines):
    logger.info('Getting xyzs and energy from XTB output')

    opt_converged, geom_section = False, False
    opt_xyzs = []
    energy = 0.0

    for line in out_lines:
        if 'GEOMETRY OPTIMIZATION CONVERGED' in line:
            opt_converged = True

        if '$coord' in line and opt_converged:
            geom_section = True

        if '$end' in line and geom_section:
            geom_section = False

        if len(line.split()) == 4 and geom_section:
            x, y, z, atom_label = line.split()
            opt_xyzs.append([atom_label,
                             float(x) * Constants.a02ang, float(y) * Constants.a02ang, float(z) * Constants.a02ang])

        if 'total E' in line:
            energy = float(line.split()[-1])

    if len(opt_xyzs) == 0:
        logger.warning('Could not find any xyzs in XTB output')
    if energy == 0:
        logger.warning('Could not find energy in XTB output')

    return opt_xyzs, energy
