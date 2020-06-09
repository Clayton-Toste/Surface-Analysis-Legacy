'''
Some real world examples for atom properties.

This is a modification of
 - psico.editing.rmsf2b
 - psico.viewing.spectrumany
to use atom properties instead of the b-factor column.
'''

import time
from pymol import cmd, CmdException

def rmsf(selection='name CA', quiet=1):
    '''
DESCRIPTION

    Determine the root mean square fluctuation (RMSF) per atom for a
    multi-state object and load it as atom property "rmsf".

ARGUMENTS

    selection = string: atom selection {default: name CA}

SEE ALSO

    spheroid, rmsf_states.py from Robert Campbell
    '''
    from numpy import array, sqrt

    n_atoms = cmd.count_atoms(selection)
    n_states = cmd.count_states(selection)
    if n_atoms == 0 or n_states < 2:
        print(' Error: not enough atoms or states')
        raise CmdException
    coords = []
    for state in range(1, n_states + 1):
        state_coords = cmd.get_model(selection, state).get_coord_list()
        if len(state_coords) != n_atoms:
            print(' Error: number of atoms in states not equal')
            raise CmdException
        coords.append(state_coords)
    coords = array(coords)
    rmsf_array = sqrt(coords.var(0).sum(1))

    cmd.alter(selection, 'properties["rmsf"] = rmsf_iter.next()',
            space={'rmsf_iter': rmsf_array.flat})

    if not int(quiet):
        print(' Average RMSF: %.2f' % (rmsf_array.mean()))
    return rmsf_array

def spectrumproperty(propname, color_list, selection='(all)', minimum=None, maximum=None, quiet=1):
    '''
DESCRIPTION

    Define a color spectrum with as many color-stops as you like (at least 2).

ARGUMENTS

    propname = string: property name which has a float value.

    color_list = string: Space separated list of colors

    ... all other arguments like with `spectrum` command
    '''
    quiet = int(quiet)
    colors = color_list.split()
    if len(colors) < 2:
        print('failed! please provide at least 2 colors')
        return

    colvec = [cmd.get_color_tuple(i) for i in colors]
    parts = len(colvec) - 1

    value_list = []
    cmd.iterate(selection, 'value_list.append(properties[propname])', space=locals())

    if len(value_list) == 0:
        print('empty selection')
        return

    if minimum is None:
        minimum = min(value_list)
    if maximum is None:
        maximum = max(value_list)
    minimum, maximum = float(minimum), float(maximum)
    val_range = (maximum - minimum) * (1 + 1e-6)

    if not quiet:
        print(' Spectrum: range (%.5f to %.5f)' % (minimum, maximum))

    if maximum == minimum:
        print('no spectrum possible, only equal values')
        return

    rgb = lambda i, p, j: int(255 * (colvec[i+1][j] * p + colvec[i][j] * (1.0 - p)))

    col_list = []
    for value in value_list:
        p = (value - minimum) / val_range * parts
        i = int(p)
        p -= i
        col_list.append(0x40000000 +
                rgb(i, p, 0) * 0x10000 +
                rgb(i, p, 1) * 0x100 +
                rgb(i, p, 2))

    cmd.alter(selection, 'color = col_iter.next()', space={'col_iter': iter(col_list)})
    cmd.recolor()

class timing(object):
    def __enter__(self):
        self.start = time.time()
    def __exit__(self, exc_type, exc_value, traceback):
        if not exc_type:
            print(' timing: %f sec' % (time.time() - self.start))

def example():
    cmd.fetch('1d7q')
    cmd.show_as('sticks')
    with timing():
        rmsf('all')
    with timing():
        spectrumproperty('rmsf', 'blue white red')

example()
