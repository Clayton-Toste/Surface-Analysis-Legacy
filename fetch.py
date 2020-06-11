from sys import argv, exit, path
from pymol import cmd, stored

'''
Script if run by main.py using the PyMOL Python Interprepeter to interface wtih PyMol.
Either loads argv[1] as protein code or file depending on argv[2]
Then colors it before saving it as 3D model to be processed
'''

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

# argv[2] determines whether argv[1] is a protein code or an file path to the protein
if (argv[2] == 'True'):
    cmd.fetch(argv[1])
else:
    cmd.load(argv[1])

#TODO
#if len(argv) == 4:
#    import pymol_filter

cmd.show_as('surface', cmd.get_object_list('all')[0])
# Set the color of each atom to match classification determined by 3-letter code of amino acid(resn)
cmd.alter('all', 'color=stored.classification_dict.get(resn, 0)')
cmd.recolor()
# Save colored 3D model of protein as .wrl file
cmd.save(cmd.get_object_list('all')[0]+".wrl", cmd.get_object_list('all')[0])