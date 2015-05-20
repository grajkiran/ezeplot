# encoding: utf-8
##############################################################################
#                                                                            #
#                  Ezeplot - Dynamical systems visualisation                 #
#                                                                            #
#   Copyright (C) 2015 Raj Kiran Grandhi <rajkiran@aero.iitkgp.ernet.in>     #
#                                                                            #
#   This program is free software: you can redistribute it and/or modify     #
#   it under the terms of the GNU General Public License as published by     #
#   the Free Software Foundation, either version 3 of the License, or        #
#   (at your option) any later version.                                      #
#                                                                            #
#   This program is distributed in the hope that it will be useful,          #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#   GNU General Public License for more details.                             #
#                                                                            #
#   You should have received a copy of the GNU General Public License        #
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                            #
##############################################################################

try:
    import tkinter as tk
    from tkinter.messagebox import showinfo
    from tkinter.font import Font
except:
    import Tkinter as tk
    from tkMessageBox import showinfo
    from tkFont import Font

import matplotlib.widgets
import dynsystem
import numpy as np
import sys, traceback
import presets
import logging
import webbrowser
import os.path

class VEntry(tk.Entry):
    # status is a static attribute that should be set to a StatusLabel widget
    # in the application code.
    status = None
    """A Validating Entry widget."""
    def __init__(self, master, textvariable, validator = None, errvariable = None,
            debug = True, command = None, **kwargs):
        """
        textvariable should be a tk.IntVar, tk.DoubleVar or tk.StringVar
        validator is the function that is called every time the content changes.
        """
        self.var = textvariable
        tk.Entry.__init__(self, master, textvariable = self.var, **kwargs)
        self.validator = validator
        self.errvariable = errvariable or tk.StringVar(self, "")
        self.debug = debug
        self.trace_id = self.var.trace('w', self.__validate)
        self.bind('<Destroy>', self.__on_destroy)
        self.default_bg = self['bg']
        if callable(command):
            self.bind('<Return>', command)
        self.set = self.var.set

    def __on_destroy(self, *args):
        # Remove the added trace when the widget is destroyed
        self.var.trace_vdelete('w', self.trace_id)

    def __validate(self, *args):
        try:
            if callable(self.validator):
                self.validator(self.get())
            else:
                self.var.get()
            self['bg'] = self.default_bg
            self.errvariable.set("")
            if self.status is not None:
                self.status.clear()
        except BaseException as e:
            self['bg'] = '#ff5555'
            self.errvariable.set(str(e))
            if self.status is not None:
                self.status.error(self.errvariable.get())
            logging.error("%s" % self.errvariable.get())
            logging.debug(traceback.format_exc())

def license_dialog(master, icon = None):
    license_text = """
Ezeplot is free software: you can
redistribute it and/or modify it under
the terms of the GNU General Public
License as published by the Free
Software Foundation, either version 3 of
the License, or (at your option) any
later version.

This program is distributed in the hope
that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the
GNU General Public License along with
this program.  If not, see
<http://www.gnu.org/licenses/>.
    """
    return showinfo("Ezeplot license", message = "Ezeplot 1.0", detail = license_text,
            parent = master)

class AboutDialog(tk.Toplevel):
    def __init__(self, master, icon = None, splash = False):
        tk.Toplevel.__init__(self, master, bg = 'black')
        if icon is not None:
            self.tk.call("wm", "iconphoto", self._w, icon)
        contents="""Ezeplot 1.0
Raj Kiran Grandhi,
Dept. of Aerospace Engineering
IIT Kharagpur - 721302
West Bengal
rajkiran@aero.iitkgp.ernet.in
http://ezeplot.example.com
"""
        try:
            icon_file = os.path.join(os.path.dirname(__file__), 'icon.ppm')
            self.icon = tk.PhotoImage(file = icon_file)
        except:
            self.icon = tk.PhotoImage()
        img = tk.Label(self, image = self.icon, bg = 'black')
        img.grid(columnspan = 2)
        frame = tk.Frame(self, bg = 'black')
        frame.grid(sticky = tk.E + tk.W, columnspan = 2)
        frame.columnconfigure(0, weight=1)
        font_b15 = Font(size = 15, weight = 'bold')
        font_b12 = Font(size = 12, weight = 'bold')
        font_reg = Font(size = 10)
        font_url = Font(size = 12, slant = 'italic', underline = 1)
        font_email = Font(size = 12, slant = 'italic')
        self.add_line(frame, text = "Ezeplot 1.0", font = font_b15)
        self.add_line(frame, text = "Raj Kiran Grandhi", font = font_b12)
        self.add_line(frame, text = "Department of Aerospace Engineering", font = font_reg)
        self.add_line(frame, text = "Indian Institute of Technology", font = font_reg)
        self.add_line(frame, text = "Kharagpur - 721302", font = font_reg)
        self.add_line(frame, text = "West Bengal, India", font = font_reg)
        self.add_line(frame, text = "rajkiran@aero.iitkgp.ernet.in", font = font_email)
        url = self.add_line(frame, text = "http://ezeplot.example.com", pady = 8,
                font = font_url, cursor = 'hand2')
        url.bind('<Button>', lambda *args: webbrowser.open(url['text']))
        #text = tk.Label(self, bg = 'black', fg = 'yellow', text = contents, justify = tk.CENTER)
        #text.grid(columnspan = 2)
        self.title('About Ezeplot')
        self.protocol('WM_DELETE_WINDOW', self.destroy)
        self.bind('<Escape>', lambda *args: self.destroy())
        img.bind('<Button>', lambda *args: self.destroy())
        self.focus_set()
        self.grab_set()
        self.transient(master)
        if splash:
            self.center_on_screen()
        else:
            btn = tk.Label(self, text = "License", font = font_b12, cursor = 'hand1',
                    bg = 'black', fg = 'white', pady = 10)
            btn.bind('<Button>', lambda arg: license_dialog(self))
            btn.grid(row = 3, column = 0)
            btn_close = tk.Label(self, text = "Close", font = font_b12, cursor = 'hand1',
                    bg = 'black', fg = 'white', pady = 10)
            btn_close.bind('<Button>', lambda *args: self.destroy())
            btn_close.grid(row = 3, column = 1)
            self.center_on_parent(master)
            self.wait_window(self)

    def add_line(self, master, text, **kwargs):
        #font = Font(**kwargs)
        label = tk.Label(master, text = text, bg = 'black', fg = 'gray90',
                **kwargs)
        label.grid(sticky = tk.W + tk.E)
        return label

    def center_on_screen(self):
        w = 250
        h = 360
        pw = self.winfo_screenwidth()
        ph = self.winfo_screenheight()
        self.geometry("%+d%+d" % ((pw-w)//2, (ph-h)//2))

    def center_on_parent(self, parent):
        w = 250
        h = 400
        pw = parent.winfo_width()
        ph = parent.winfo_height()
        px = parent.winfo_rootx()
        py = parent.winfo_rooty()
        x = px + (pw-w)//2
        y = py + (ph-h)//2
        geom = "%+d%+d" % (x,y)
        self.geometry(geom)

class PlotLimits(tk.Toplevel):
    def __init__(self, master, fig, limits):
        tk.Toplevel.__init__(self, master)
        self.fig = fig
        self.limits = limits
        frame = tk.Frame(self)
        frame.pack()
        self.title("Change limits")
        self.transient(master)
        #tk.Label(frame, text = "Axis").grid(row = 0, column = 0)
        tk.Label(frame, text = "Min").grid(row = 0, column = 1)
        tk.Label(frame, text = "Max").grid(row = 0, column = 2)
        tk.Label(frame, text = "x:").grid(row = 1, column = 0)
        VEntry(frame, textvariable = limits.xmin, width = 5).grid(row = 1, column = 1)
        VEntry(frame, textvariable = limits.xmax, width = 5).grid(row = 1, column = 2)
        tk.Label(frame, text = "y:").grid(row = 2, column = 0)
        VEntry(frame, textvariable = limits.ymin, width = 5).grid(row = 2, column = 1)
        VEntry(frame, textvariable = limits.ymax, width = 5).grid(row = 2, column = 2)
        tk.Label(frame, text = "z:").grid(row = 3, column = 0)
        VEntry(frame, textvariable = limits.zmin, width = 5).grid(row = 3, column = 1)
        VEntry(frame, textvariable = limits.zmax, width = 5).grid(row = 3, column = 2)
        tk.Label(frame, text = "Computation factor:").grid(row = 4, column = 0, columnspan = 2, sticky = tk.E)
        VEntry(frame, textvariable = limits.factor, width = 5).grid(row = 4, column = 2)
        #tk.Label(frame, text = "Periodicty:").grid(row = 4, columnspan = 3, sticky = tk.W)
        #tk.Checkbutton(frame, text = "x", variable = limits.per_x).grid(row = 5, column = 0)
        #tk.Checkbutton(frame, text = "y", variable = limits.per_y).grid(row = 5, column = 1)
        #tk.Checkbutton(frame, text = "z", variable = limits.per_z).grid(row = 5, column = 2)
        tk.Button(frame, text = "From plot", command = self.__from_plot).grid(row = 6, column = 0, columnspan = 2)
        tk.Button(frame, text = "Close", command = self.destroy).grid(row = 6, column = 2)
        self.protocol('WM_DELETE_WINDOW', self.destroy)
        self.grab_set()
        self.wait_window(self)

    def __from_plot(self, *args):
        limits = self.fig.get_limits()
        xmin, xmax = limits[0]
        ymin, ymax = limits[1]
        zmin, zmax = limits[2]
        self.limits.xmin.set(xmin)
        self.limits.xmax.set(xmax)
        self.limits.ymin.set(ymin)
        self.limits.ymax.set(ymax)
        self.limits.zmin.set(zmin)
        self.limits.zmax.set(zmax)

class PEntry(tk.Frame):
    def __init__(self, master, name, textvariable = None, validator = None,
            command = None, value = 1.0):
        tk.Frame.__init__(self, master)
        self.name = name
        self.var = textvariable or tk.DoubleVar(self, value)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.label = tk.Label(self, text = name, width = 6)
        self.entry = VEntry(self, self.var, width = 6, validator = validator,
                command = command)
        self.get = self.var.get
        self.set = self.var.set
        self.enable()
    def enable(self, name = None, value = None):
        if name is not None:
            self.name = name
        self.label.configure(text = self.name)
        if value is not None:
            self.var.set(value)
        self.label.grid(sticky = tk.E + tk.W)
        self.entry.grid(row = 0, column = 1, sticky = tk.E+tk.W)
    def disable(self):
        self.label.grid_forget()
        self.entry.grid_forget()

class DSFrame(tk.LabelFrame):
    def __init__(self, master, system, command = None, n_params = 10, preset_cmd = None, **kwargs):
        tk.LabelFrame.__init__(self, master, text = "Dynamical system", **kwargs)
        self.columnconfigure(1, weight=1)
        self.system = system
        self.command = command
        self.preset_cmd = preset_cmd
        self.eqn_x = tk.StringVar(self, self.system.get(1))
        self.eqn_y = tk.StringVar(self, self.system.get(2))
        self.eqn_z = tk.StringVar(self, self.system.get(3))
        self.preset = tk.StringVar(self)
        choices = list(presets.systems.keys())
        #choices.append("User defined")
        tk.OptionMenu(self, self.preset, *choices,
                command = self._load_preset).grid(sticky = tk.E+ tk.W, columnspan = 2)
        tk.Label(self, text = "x_dot:").grid()
        self.entry_x = VEntry(self, self.eqn_x, command = self.command,
                validator = self._update_system, disabledforeground = "black")
        self.entry_x.grid(row=1, column = 1, sticky = tk.E + tk.W)
        tk.Label(self, text = "y_dot:").grid()
        self.entry_y = VEntry(self, self.eqn_y, command = self.command,
                validator = self._update_system, disabledforeground = "black")
        self.entry_y.grid(row=2, column = 1, sticky = tk.E + tk.W)
        tk.Label(self, text = "z_dot:").grid()
        self.entry_z = VEntry(self, self.eqn_z, command = self.command,
                validator = self._update_system, disabledforeground = "black")
        self.entry_z.grid(row=3, column = 1, sticky = tk.E + tk.W)
        self.params = list()
        pframe = tk.Frame(self)
        pframe.grid(columnspan = 2, sticky = tk.E + tk.W)
        pframe.columnconfigure(0, weight=1)
        pframe.columnconfigure(1, weight=1)
        for row in range(n_params//2):
            for col in 0, 1:
                pe = PEntry(pframe, "", value = 0.0, command = self.command,
                        validator = self._update_system_params)
                pe.grid(row = row, column = col, sticky = tk.E + tk.W)
                self.params.append(pe)
        #self._load_preset(choices[-1])

    def _load_preset(self, name):
        self.preset.set(name)
        preset = presets.systems[name]
        #x_dot, y_dot, z_dot, params = dynsystem.presets[name]
        x_dot = preset['x']
        y_dot = preset['y']
        z_dot = preset['z']
        params = preset.get('params', dict())
        self.eqn_x.set(x_dot)
        self.eqn_y.set(y_dot)
        self.eqn_z.set(z_dot)
        self.system.params.update(params)
        self._update_params()
        if callable(self.preset_cmd):
            self.preset_cmd(name)
        if name != "User defined":
            for e in self.entry_x, self.entry_y, self.entry_z:
                e.configure(state = tk.DISABLED)
        else:
            for e in self.entry_x, self.entry_y, self.entry_z:
                e.configure(state = tk.NORMAL)
        self._update_params()

    def _update_system(self, *args):
        x_dot = self.eqn_x.get()
        y_dot = self.eqn_y.get()
        z_dot = self.eqn_z.get()
        self.system.update(x_dot, y_dot, z_dot)
        self._update_params()

    def _update_system_params(self, *args):
        for pe in self.params:
            name = pe.name
            value = pe.get()
            if name in self.system.params:
                self.system.params[name] = value

    def _update_params(self, *args):
        items = list(self.system.params.items())
        for i in range(len(items)):
            self.params[i].enable(name = items[i][0], value = items[i][1])
        for pe in self.params[len(items):]:
            pe.name = ""
            pe.disable()

class StatusLabel(tk.Label):
    def __set_message(self, message):
        self.configure(text = message)
        self.update_idletasks()

    def error(self, message):
        self['fg'] = "#aa0000"
        self.__set_message(message)

    def warn(self, message):
        self['fg'] = "#ff5555"
        self.__set_message(message)

    def info(self, message):
        self['fg'] = "black"
        self.__set_message(message)

    def ready(self):
        self['fg'] = "black"
        self.__set_message("Ready.")

    def clear(self):
        self['fg'] = 'black'
        self.__set_message("")

def dsf_test():
    root = tk.Tk()
    d = dynsystem.DynamicSystem("y+a*y-b*x", "mu*x-c*y", dict(mu = 1.0))
    f = DSFrame(root, d)
    root.columnconfigure(0, weight=1)
    f.grid(sticky = tk.E+tk.W, columnspan=3)
    e1 = FloatEntry(root, default = 1.0)
    e2 = FloatEntry(root, default = 2.0)
    b = tk.Button(root, command = test_system, text = "Try")
    e1.grid(row=1, column = 0)
    e2.grid(row=1, column = 1)
    b.grid(row=1, column = 2)
    root.mainloop()
if __name__ == '__main__':
    root = tk.Tk()
    ad = AboutDialog(root)
    #import IPython
    #IPython.embed()
