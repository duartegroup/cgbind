"""

Building blocks to generate a library of linkers suitable to form [M2L4]n+ cages, with a general form


        ------top     --- M
        |
        |
        link
        |
        |
        centre
        |
        |
        link
        |
        |
        ------bottom  --- M

"""


def make_linker_smiles(top_smiles, top_link_smiles, centre_smiles, bottom_link_smiles, bottom_smiles):
    return top_smiles + top_link_smiles + centre_smiles + bottom_link_smiles + bottom_smiles


top_end_names_smiles = {
    '3-pyridine'     : 'C1=CC=NC=C1%99',
    '7-isoquinoline' : 'C1%99=CC2=C(C=C1)C=CN=C2',

}

bottom_end_names_smiles = {
    '3-pyridine'     : 'C1%96=CC=CN=C1',
    '7-isoquinoline' : 'C1%96=CC2=C(C=C1)C=CN=C2'
}

top_link_names_smiles = {
    'alkyne'        : '.C%99#C%98.',
    '1-4-phenyl'    : '.C1%99=CC=C%98C=C1.',
    '1-4-napythyl'  : '.C1%99=C(C=CC=C2)C2=C%98C=C1.',
    '2-4-BCP'       : '.C1%99(C2)CC2%98C1.',
    'NH'            : '.N%99%98.',
    'amide1'        : '.C%99(=O)N%98.',
    'amide2'        : '.N%99C%98(=O).',
    'phenyl_alkyne1': '.C1%99=CC=C(C#C%98)C=C1.',
    'phenyl_alkyne2': '.C%99#CC1=CC=C%98C=C1.',
    'piperazine'    : '.N%991CCN%98CC1.'
}

bottom_link_names_smiles = {
    'alkyne'        : '.C%97#C%96.',
    '1-4-phenyl'    : '.C1%97=CC=C%96C=C1.',
    '1-4-napythyl'  : '.C1%97=C(C=CC=C2)C2=C%96C=C1.',
    '2-4-BCP'       : '.C1%97(C2)CC2%96C1.',
    'NH'            : '.N%97%96.',
    'amide1'        : '.C%97(=O)N%96.',
    'amide2'        : '.N%97C%96(=O).',
    'phenyl_alkyne1': '.C1%97=CC=C(C#C%96)C=C1.',
    'phenyl_alkyne2': '.C%97#CC1=CC=C%96C=C1.',
    'piperazine'    : '.N%971CCN%96CC1.'
    }

centre_names_smiles = {
    '2-6-pyridine'         : 'C1%98=NC%97=CC=C1',
    'carbonate'            : 'O=C(O%98)O%97',
    '2-7-napthaly'         : 'C1%98=CC=C(C=CC%97=C2)C2=C1',
    '1-8-anthracene'       : 'C1%98=CC=CC2=C1C=C3C(C=CC=C3%97)=C2',
    '1-3-phenyl'           : 'C1%98=CC=CC%97=C1',
    '1-3-phenyl-2-methyl'  : 'C1%98=CC=CC%97=C1C',
    'methylene'            : 'C%98%97',
}