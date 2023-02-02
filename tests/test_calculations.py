from cgbind.calculations import _cgbind_mol_to_autode
from cgbind.molecule import BaseStruct, Molecule
import numpy as np


def check_equivalence_ade_and_cgbind_mol(ade_mol, cgb_mol):
    assert cgb_mol.charge == ade_mol.charge
    assert cgb_mol.mult == ade_mol.mult
    if cgb_mol.solvent is not None:
        assert cgb_mol.solvent == str(ade_mol.solvent)
    else:
        assert ade_mol.solvent is None
    assert cgb_mol.name == ade_mol.name
    assert cgb_mol.n_atoms == ade_mol.n_atoms

    coords1 = np.array(ade_mol.coordinates)
    coords2 = np.array(cgb_mol.get_coords())

    assert np.isclose(coords1, coords2).all()


def test_cgbind_to_ade_mol_conversion():
    mol = Molecule(smiles='CC[O-]', name='test', charge=-1,
                   mult=1, solvent='water')

    ade_mol = _cgbind_mol_to_autode(mol)
    check_equivalence_ade_and_cgbind_mol(ade_mol, mol)


def test_cgbind_to_ade_conformer_conversion():
    # only works with SMILES as of now (2 Feb 2023)
    mol = Molecule(smiles='CCCCCO', name='test', n_confs=5)

    # n_confs is not the true number!!
    ade_mol = _cgbind_mol_to_autode(mol)
    assert ade_mol.n_conformers == len(mol.conformers) - 1

    for idx in range(ade_mol.n_conformers):
        cgb_conf = mol.conformers[idx+1]
        ade_conf = ade_mol.conformers[idx]
        check_equivalence_ade_and_cgbind_mol(ade_conf, cgb_conf)
