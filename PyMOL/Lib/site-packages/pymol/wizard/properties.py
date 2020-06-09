from pymol import wizard, querying

class Properties(wizard.Wizard):

    obj_name = None
    atom_idx = 0
    prop_list = None

    default_prompt = ['Please click an object...']
    prompt = default_prompt

    mode = 0
    mode_name = [
        'Object',
        'Atom',
    ]

    prop = 'all'

    def __init__(self, *args, **kwargs):
        wizard.Wizard.__init__(self, *args, **kwargs)

        self.cmd.set("mouse_selection_mode", 0)

        self.menu = {
            'mode': [
                [2, 'Level', ''],
            ] + [
                [1, label, 'cmd.get_wizard().set_mode(%d)' % i]
                for (i, label) in enumerate(self.mode_name)
            ],
            'prop': [
                [2, 'Property', ''],
                [1, 'all', 'cmd.get_wizard().set_prop("all")'],
                [0, '', ''],
            ],
        }

    def get_menu(self, key):
        if key == 'prop':
            self.menu[key][3:] = [
                [1, label, 'cmd.get_wizard().set_prop(%s)' % repr(label)]
                for label in (self.prop_list or [])
            ]
        return self.menu[key]

    def get_event_mask(self):
        return \
            self.event_mask_pick + \
            self.event_mask_select + \
            self.event_mask_state

    def get_panel(self):
        return [
            [ 1, 'Properties',''],
            [ 3, 'Level: ' + self.mode_name[self.mode], 'mode'],
            [ 3, 'Property: ' + self.prop, 'prop'],
            [ 2, 'Clear','cmd.get_wizard().clear_prompt()'],
            [ 2, 'Done','cmd.set_wizard()'],
            ]

    def set_mode(self, mode):
        self.mode = mode
        self.clear_prompt()

    def set_prop(self, prop):
        self.prop = prop
        self.cmd.refresh_wizard()

    def clear_prompt(self):
        self.obj_name = None
        self.prop = 'all'
        self.cmd.refresh_wizard()

    def do_state(self, state):
        self.cmd.refresh_wizard()

    def do_select(self, name):
        self.cmd.edit('first (%s)' % (name))
        self.do_pick()

    def do_pick(self, bondFlag=0):
        if self.cmd.count_atoms('?pk2'):
            self.cmd.select('pk1', 'pk2')
            self.cmd.delete('pk2')
        self.obj_name, self.atom_idx = self.cmd.index('pk1')[0]
        self.cmd.refresh_wizard()

    def get_prompt(self):
        self.update_prompt()
        return self.prompt

    def _get_property(self, key, name, state):
        if self.mode:
            return self.prop_list.get(key)
        return self.cmd.get_property(key, name, state)

    def update_prompt(self):
        if not self.obj_name:
            if self.mode:
                self.prompt = ['Please click an atom...']
            else:
                self.prompt = self.default_prompt
            return

        name = self.obj_name
        state = querying.get_object_state(name)

        self.prompt = [
            r'\990Object: %s, State: %d' % (name, state),
            r'\777',
        ]

        if self.mode:
            sele = '%s`%d' % (name, self.atom_idx)
            n = self.cmd.iterate_state(state, sele, 'self.atom_macro = "/%s/%s/%s/%s`%s/%s" % '
                    '(model, segi, chain, resn, resi, name);'
                    'self.prop_list = properties.all', space={'self': self})
            if n == 0:
                raise ValueError
            self.prompt[0] = r'\990Atom: %s, State: %d' % (self.atom_macro, state)
        else:
            self.prop_list = self.cmd.get_property_list(name, state)

        if self.prop != 'all':
            value = self._get_property(self.prop, name, state)
            self.prompt.append(r'\669%s:\777' % (self.prop))
            if isinstance(value, (int, float)):
                self.prompt[-1] += ' ' + str(value)
            else:
                self.prompt.append('')
                self.prompt.extend(str(value).splitlines())
            return

        if not self.prop_list:
            self.prompt += [
                r'\955- empty property list -',
            ]
            return

        tlen = 20 + max(map(len, self.prop_list))

        for i, key in enumerate(self.prop_list):
            value = self._get_property(key, name, state)

            try:
                alen = tlen - len(key)
                svalue = str(value)
                if len(svalue) > alen:
                    svalue = svalue[:alen - 3] + '...'
            except:
                svalue = '?'
            self.prompt.append(r'\669%s:\777 %s' % (key, svalue))
