'''
Setting help

Reads settings documentation from a CSV file into a dictionary.
'''

from pymol import cmd

def setting_help_read_csv(filename):
    import csv
    csvreader = csv.reader(open(filename), escapechar='\\', doublequote=False)
    return dict((row[0], row[1:]) for row in csvreader)

setting_help_dict = setting_help_read_csv(
        cmd.exp_path('$PYMOL_DATA/setting_help.csv'))
