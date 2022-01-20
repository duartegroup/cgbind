from cgbind.utils import fast_xtb_opt
from cgbind.substrate import Substrate
from cgbind.atoms import Atom
import shutil


def test_simple_opt():

    methane = Substrate(name='methane')
    methane.set_atoms(atoms=[Atom('C', -0.0221, 0.0032, 0.0165),
                             Atom('H', -0.6690, 0.8894, -0.100),
                             Atom('H', -0.3778, -0.857, -0.588),
                             Atom('H', 0.09640, -0.315, 1.0638),
                             Atom('H', 0.97250, 0.2803, -0.391)])

    fast_xtb_opt(methane, n_cores=1, n_cycles=5)

    assert methane.atoms is not None

    if shutil.which('xtb') is not None:

        # Should have shifted the atom a little
        assert methane.atoms[0].coord[0] != -0.0221
