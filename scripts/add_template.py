"""
Add a template to the library of available templates e.g.

python add_template.py metallocage_strucutre.mol2 --name name_of_this_template
"""
import argparse
from cgbind.templates import Template


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("filename", action='store', type=str,
                        help='.mol2 file) add to the template library')

    parser.add_argument("-n", '--name', action='store', type=str,
                        help='Name of the template')

    return parser.parse_args()


if __name__ == '__main__':

    args = get_args()
    template = Template(arch_name=args.name, mol2_filename=args.filename)
    template.save_template()
