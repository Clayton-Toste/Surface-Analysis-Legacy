from sys import argv, exit, path
from pymol import cmd, stored

stored.classification_dict = {
    'ASP': 19,
    'GLU': 19,
    'ARG': 3,
    'LYS': 3,
    'HIS': 3,
    'ASN': 2,
    'SER': 2,
    'GLN': 2,
    'PRO': 2,
    'THR': 2, 
    'TYR': 2, 
    'CYS': 2,
    'TRP': 4,
    'ALA': 4,
    'ILE': 4,
    'LEU': 4,
    'MET': 4, 
    'PHE': 4,
    'VAL': 4,
    'GLY': 4,
}

if (argv[2] == 'True'):
    cmd.fetch(argv[1])
else:
    cmd.load(argv[1])
if len(argv) == 4:
    import pymol_filter
cmd.hide('everything')
cmd.show('surface', cmd.get_object_list('all')[0])
cmd.alter('all', 'color=stored.classification_dict.get(resn, 0)')
cmd.recolor()
cmd.save(cmd.get_object_list('all')[0]+".wrl", cmd.get_object_list('all')[0])