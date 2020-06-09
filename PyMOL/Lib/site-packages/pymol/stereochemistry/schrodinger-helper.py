##
## Schrodinger helper script, prints R/S stereo labels, one line per atom
##

if __name__ != '__main__':
    raise ImportError('run with $SCHRODINGER/run')

import sys
from schrodinger import structure
from schrodinger.structutils import analyze

for filename in sys.argv[1:]:
    for st in structure.StructureReader(filename):
        chiralities = analyze.get_chiral_atoms(st)
        for a in st.atom:
            print(chiralities.get(a.index, ''))
