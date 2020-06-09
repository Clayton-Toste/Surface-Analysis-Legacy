from pymol.wizard.command import *

class Morphing(Command):

    ignored_args = ('quiet', 'name', 'match')

    def __init__(self, sele1='', _self=cmd):
        Command.__init__(self, 'morph', _self)

        self.current['name'] = ''
        self.current['sele1'] = sele1
        self.current['sele2'] = ''

        self.set_menu_values('method', ('rigimol', 'linear'), last=None)
        self.set_menu_values('refinement', [0, 1, 2, 5, 10])
        self.set_menu_values('steps', [10, 20, 30, 50, 70, 100])

        for arg in ('state1', 'state2'):
            self.menu[arg][1:] = [
                [1, '-1 (auto)', '%s.set_current("%s", -1)' % (self.varname, arg) ],
                [1, '0 (all)' if arg == 'state1' else '0 (last)',
                    '%s.set_current("%s",  0)' % (self.varname, arg) ],
            ]

    def get_menu(self, arg):
        if arg in ('state1', 'state2'):
            try:
                value = self.current['sele' + arg[-1]]
                nstates = self.cmd.count_states(value)
            except:
                nstates = 0
            self.set_menu_values(arg, (range(1, 9) + [nstates] if nstates > 10
                else range(1, nstates+1)), 3, None)
        return Command.get_menu(self, arg)
