import pymol

class Atomprops(pymol.wizard.Wizard):

    temp_settings = [
        ("mouse_selection_mode", 0),
        ("overlay_lines", 100),
        ("overlay", 1),
    ]

    def __init__(self, sele='', *args, **kwargs):
        pymol.wizard.Wizard.__init__(self, *args, **kwargs)

        self.incl_obj = 1
        self.state = -1
        self.space = {'callback_state': self.callback_state, 'callback': self.callback}
        self.orig_settings = []

        for name, value in self.temp_settings:
            current = self.cmd.get_setting_int(name)
            if value != current:
                self.orig_settings.append((name, current))
                self.cmd.set(name, value)

        self.cmd.cls()

        if sele:
            self.cmd.edit('first (' + sele + ')')
            self.do_pick()

    def cleanup(self):
        for name, value in self.orig_settings:
            self.cmd.set(name, value)

    def safe_repr(self, s):
        s = repr(s)
        if len(s) > 80:
            s = s[:77] + '...'
        return s

    def get_event_mask(self):
        return \
            self.event_mask_pick + \
            self.event_mask_select + \
            self.event_mask_state

    def get_panel(self):
        return [
            [ 1, 'Atom Properties and Settings', ''],
            [ 3, 'w/ object-level: ' + ('Yes' if self.incl_obj else 'No'), 'incl_obj'],
            [ 2, 'Done','cmd.set_wizard()'],
            ]

    def get_menu(self, key):
        return [
            [1, 'Yes', 'cmd.get_wizard().set_incl_obj(True)'],
            [1, 'No', 'cmd.get_wizard().set_incl_obj(False)'],
        ]

    def set_incl_obj(self, value):
        self.incl_obj = value
        self.cmd.refresh_wizard()

    def get_prompt(self):
        return ['Please click an atom...']

    def do_state(self, state):
        self.cmd.refresh_wizard()

    def do_select(self, name):
        self.cmd.edit('first ?' + name)
        self.do_pick()

    def do_pick(self, bondFlag=0):
        cmd = self.cmd

        if cmd.count_atoms('?pk2'):
            cmd.select('pk1', '?pk2')
            cmd.delete('pk2')

        cmd.iterate('?pk1', 'callback(model, segi, chain, resn, resi, name, alt,'
                'p, s, index, rank, ID, b, q, formal_charge, partial_charge, elem,'
                'protons, type, vdw, ss, color, reps, custom, elec_radius)',
                space=self.space)

        cmd.refresh_wizard()

    def callback_state(self, s, state):
        keys = list(s)
        if keys:
            print('  settings (atom-state: %d):' % (state))
            for i in keys:
                print('    %3d %s: %s' % (i, pymol.setting.name_dict.get(i, '<unknown>'), self.safe_repr(s[i])))

    def callback(self, model, segi, chain, resn, resi, name, alt, p, s,
            index, rank, ID, b, q, formal_charge, partial_charge,
            elem, protons, type, vdw, ss, color, reps, custom, elec_radius):
        self.cmd.cls()

        color = ('%d' if color < 0x40000000 else '0x%x') % color

        print('/%s/%s/%s/%s`%s/%s`%s' % (model, segi, chain, resn, resi, name, alt))
        print('  index:%d rank:%d ID:%d b:%.2f q:%.2f vdw:%.2f' % (index, rank, ID, b, q, vdw))
        print('  formal_charge:%d partial_charge:%.2f elec_radius:%.2f' % (
            formal_charge, partial_charge, elec_radius))
        print('  elem:%s protons:%d type:%s ss:%s custom:%s color:%s reps:%s' % (
            elem, protons, type, ss, custom, color, bin(reps)))

        keys = list(p)
        if keys:
            print('  properties (atom-level):')
            for key in keys:
                print('    p.' + str(key) + ': ' + self.safe_repr(p[key]))

        keys = list(s)
        if keys:
            print('  settings (atom-level):')
            for i in keys:
                print('    %3d %s: %s' % (i, pymol.setting.name_dict.get(i, '<unknown>'), self.safe_repr(s[i])))

        self.cmd.iterate_state(self.state, (model, index),
                'callback_state(s, state)',
                space=self.space)

        if self.incl_obj:
            keys = self.cmd.get_property_list(model, self.state)
            if keys:
                print('  properties (object-state-level):')
                for key in keys:
                    value = self.cmd.get_property(key, model, self.state)
                    print('    ' + str(key) + ': ' + self.safe_repr(value))

            for state, label in [(0, 'object'), (self.state, 'object-state')]:
                s_list = self.cmd.get_object_settings(model, state)
                if s_list:
                    print('  settings (%s-level):' % label)
                    for sitem in s_list:
                        key = pymol.setting.name_dict.get(sitem[0], sitem[0])
                        print('    ' + str(key) + ': ' + self.safe_repr(sitem[2]))
