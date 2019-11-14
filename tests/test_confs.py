from rdkit.Chem import AllChem
from rdkit import Chem
from cgbind.confomers import extract_xyzs_from_rdkit_mol_object


def test_conf_conv():

    propane = Chem.MolFromSmiles('CCC')
    propane = Chem.AddHs(propane)
    conf_ids = list(AllChem.EmbedMultipleConfs(propane, numConfs=2, params=AllChem.ETKDG()))

    # Generate the .mol files

    conf_xyzs = extract_xyzs_from_rdkit_mol_object(mol_obj=propane, conf_ids=conf_ids)
    assert len(conf_xyzs) == 2
    assert len(conf_xyzs[0]) == 11      # propane has 11 atoms

