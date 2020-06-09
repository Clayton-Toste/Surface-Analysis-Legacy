'''
Callout/Label settings panel

(c) 2013 Schrodinger Inc.
'''

from __future__ import absolute_import
from __future__ import division

import sys

if sys.version_info[0] == 2:
    import Tkinter
    import tkMessageBox
else:
    import tkinter as Tkinter
    from tkinter import messagebox as tkMessageBox

import Pmw
from pymol import cmd, plugins
from pmg_tk import Setting

def sliderentry(master, var, **kwargs):
    frame = Tkinter.Frame(master)
    scale = Tkinter.Scale(frame, showvalue=0, orient=Tkinter.HORIZONTAL,
            variable=var, **kwargs)
    entry = Tkinter.Entry(frame, textvariable=var, width=5)
    entry.pack(side=Tkinter.RIGHT)
    scale.pack(side=Tkinter.LEFT, expand=1, fill=Tkinter.X)
    return frame

def entrynf(master, var, n=3, width=0, **kwargs):
    if not width:
        width = 16 // n
    frame = Tkinter.Frame(master)
    for i in range(n):
        Tkinter.Entry(frame, textvariable=var[i],
                width=width, **kwargs).pack(side=Tkinter.LEFT)
    return frame

def unset_label_settings(selection):
    # skipped: label_screen_point, label_position
    names = ["label_anchor", "label_angle_digits", "label_bg_color",
            "label_bg_outline", "label_bg_transparency", "label_color",
            "label_connector", "label_connector_color",
            "label_connector_ext_length", "label_connector_mode",
            "label_connector_width", "label_digits", "label_dihedral_digits",
            "label_distance_digits", "label_font_id",
            "label_multiline_justification", "label_multiline_spacing",
            "label_outline_color", "label_padding", "label_placement_offset",
            "label_relative_mode", "label_shadow_mode", "label_size"]
    for name in names:
        cmd.unset(name, selection)

def __init_plugin__(self=None):
    plugins.addmenuitem('Label Settings', calloutpanel)

def addrow(*widgets, **kwargs):
    for j, w in enumerate(widgets):
        w.grid(row=addrow.i, column=j, sticky=Tkinter.W, **kwargs)
    addrow.i += 1
addrow.i = 0

@cmd.extend
def calloutpanel():
    self = Tkinter.Toplevel(plugins.get_tk_root())
    self.title('Label Settings')
    master = self
    app = plugins.get_pmgapp()

    def makeframe(text):
        f = Tkinter.LabelFrame(self, text=text, padx=5, pady=5)
        f.pack(fill=Tkinter.X, expand=1)
        return f

    topframes = {}
    def updatetopframe(v):
        for top in topframes.values():
            top.pack_forget()
        if v not in topframes:
            if v == 'PICKED':
                if cmd.count_atoms('?pk1') != 1:
                    tkMessageBox.showwarning('Warning', 'Please pick an atom in edit mode')
                    return
                setting = Setting.AtomSetting(app, 'pk1')
            elif v == 'GLOBAL':
                setting = app.skin.setting
            else:
                raise ValueError(v)
            topframes[v] = maketopframe(self, setting)
        topframes[v].pack(fill=Tkinter.BOTH, expand=1)

    but_level = Pmw.OptionMenu(master, labelpos='w', label_text='Selection:',
            items=['GLOBAL', 'PICKED'], command=updatetopframe)
    but_level.pack(side=Tkinter.TOP)

    master = makeframe("Unset Object or Atom level Settings")
    Tkinter.Button(master, text="all atoms",
            command=lambda: unset_label_settings('(*)')
            ).pack(side=Tkinter.LEFT)
    Tkinter.Button(master, text="all objects",
            command=lambda: unset_label_settings('*')
            ).pack(side=Tkinter.LEFT)

    updatetopframe('GLOBAL')

def maketopframe(self, setting):
    top = Tkinter.Frame(self)

    def makeframe(text):
        f = Tkinter.LabelFrame(top, text=text, padx=5, pady=5)
        f.pack(fill=Tkinter.X, expand=1)
        return f

    ## GENERAL
    master = makeframe("General")

    # label_relative_mode           : 0
    addrow(Tkinter.Label(master, text="Relative Mode"),
           Tkinter.OptionMenu(master, setting.label_relative_mode, "0", "1", "2"))

    # float_labels : 0
    addrow(Tkinter.Label(master, text="Floating Labels"),
           Tkinter.Checkbutton(master, variable=setting.float_labels))

    # label_position                : [ 0.00000, 0.00000, 1.75000 ]
    addrow(Tkinter.Label(master, text="Position"),
            entrynf(master, setting.label_position))

    ### CONNECTOR
    master = makeframe("Connector")

    # label_connector               : off
    but_connector = Tkinter.Checkbutton(master, variable=setting.label_connector)

    # label_connector_mode          : 0
    but_connector_mode = sliderentry(master, setting.label_connector_mode,
            from_=0, to=4, resolution=1)

    # label_connector_color         : front
    but_connector_color = Tkinter.OptionMenu(master, setting.label_connector_color,
            "front", "white", "gray", "black")

    # label_connector_width         : 2.00000
    but_connector_width = sliderentry(master, setting.label_connector_width,
            from_=0.0, to=8.0, resolution=0.1)

    # label_connector_ext_length    : 2.50000
    but_connector_ext_length = sliderentry(master, setting.label_connector_ext_length,
            from_=0.0, to=8.0, resolution=0.1)

    addrow(Tkinter.Label(master, text="Show Connector"), but_connector)
    addrow(Tkinter.Label(master, text="Mode"), but_connector_mode)
    addrow(Tkinter.Label(master, text="Color"), but_connector_color)
    addrow(Tkinter.Label(master, text="Line Width"), but_connector_width)
    addrow(Tkinter.Label(master, text="Extension Length"), but_connector_ext_length)

    ### BACKGROUND
    master = makeframe("Background")

    # label_bg_color                : default
    addrow(Tkinter.Label(master, text="Color"),
            Tkinter.OptionMenu(master, setting.label_bg_color,
                "default", "back", "front", "white", "gray", "black"))

    # label_bg_transparency         : 0.60000
    addrow(Tkinter.Label(master, text="Transparency"),
            sliderentry(master, setting.label_bg_transparency,
                from_=0.0, to=1.0, resolution=0.1))

    # label_bg_outline              : off
    addrow(Tkinter.Label(master, text="Show Outline"),
            Tkinter.Checkbutton(master, variable=setting.label_bg_outline))

    # label_padding                 : [ 0.20000, 0.20000, 0.00000 ]
    addrow(Tkinter.Label(master, text="Padding"),
            entrynf(master, setting.label_padding, 2))

    ### FONT
    master = makeframe("Font")

    # label_color                   : front
    but_label_color = Tkinter.OptionMenu(master, setting.label_color,
            "front", "back", "atomic", "white", "gray", "black")

    # label_size                    : 14.00000
    but_label_size = sliderentry(master, setting.label_size, from_=5.0, to=50.0, resolution=0.1)

    addrow(Tkinter.Label(master, text="Color"), but_label_color)
    addrow(Tkinter.Label(master, text="Font Size"), but_label_size)

    # label_multiline_justification : 1.00000
    addrow(Tkinter.Label(master, text="Justification"),
            sliderentry(master, setting.label_multiline_justification,
                from_=1, to=-1, resolution=1))
                #from_=-1, to=1, resolution=1))

    # label_multiline_spacing       : 1.20000
    addrow(Tkinter.Label(master, text="Line Spacing"),
            sliderentry(master, setting.label_multiline_spacing,
                from_=1.0, to=2.0, resolution=0.1))
    

    # label_font_id                 : 5
    var_font_family = Tkinter.StringVar(master, "Sans")
    var_font_bold = Tkinter.IntVar(master)
    var_font_italic = Tkinter.IntVar(master)

    def update_font_id(*args, **kwargs):
        family = var_font_family.get()
        style = 1 + var_font_bold.get() + var_font_italic.get() * 2
        if family == 'Serif':
            font_id = 10 if style == 2 else \
                      17 if style == 3 else \
                      18 if style == 4 else \
                      9
        elif family == 'Mono':
            font_id = 13 if style == 2 else \
                      12 if style == 3 else \
                      14 if style == 4 else \
                      11
        else:
            font_id = 7 if style == 2 else \
                      6 if style == 3 else \
                      8 if style == 4 else \
                      5
        setting.label_font_id.set(font_id)

    var_font_family.trace("w", update_font_id)
    var_font_bold.trace("w", update_font_id)
    var_font_italic.trace("w", update_font_id)

    font_id_frame = Tkinter.Frame(master)
    Tkinter.OptionMenu(font_id_frame, var_font_family,
            "Sans", "Serif", "Mono").pack(side=Tkinter.LEFT)
    Tkinter.Checkbutton(font_id_frame, text="Bold", variable=var_font_bold).pack(
            side=Tkinter.LEFT)
    Tkinter.Checkbutton(font_id_frame, text="Italic", variable=var_font_italic).pack(
            side=Tkinter.LEFT)
    addrow(Tkinter.Label(master, text="Family"), font_id_frame)

    return top

# vi:expandtab:smarttab
