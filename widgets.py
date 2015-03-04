#!/usr/bin/python3

try:
    import tkinter as tk
except:
    import Tkinter as tk

import matplotlib.widgets
import dynsystem
import numpy as np
from traceback import print_exc
import sys

class VEntry(tk.Entry):
    """A Validating Entry widget."""
    def __init__(self, master, textvariable, validator = None,
            debug = True, command = None, **kwargs):
        """
        textvariable should be a tk.IntVar, tk.DoubleVar or tk.StringVar
        validator is the function that is called every time the content changes.
        """
        self.var = textvariable
        tk.Entry.__init__(self, master, textvariable = self.var, **kwargs)
        self.validator = validator
        self.errmsg = tk.StringVar(self, "")
        self.debug = debug
        self.var.trace('w', self.__validate)
        self.default_bg = self['bg']
        if callable(command):
            self.bind('<Return>', command)
        self.set = self.var.set

    def __validate(self, *args):
        try:
            if callable(self.validator):
                self.validator(self.get())
            else:
                self.var.get()
            self['bg'] = self.default_bg
            self.errmsg.set("")
        except BaseException as e:
            self['bg'] = '#ff5555'
            self.errmsg.set(str(e))
            if self.debug:
                sys.stderr.write("%s\n" % self.errmsg.get())

class Param(tk.LabelFrame):
    def __init__(self, master, name, value = 1.0, command = None):
        tk.LabelFrame.__init__(self, master, text = name)
        self.name = name
        self.command = command
        self.columnconfigure(0, weight=1)
        self.entry = VEntry(self, tk.StringVar(master, str(value)), debug = True,
                validator = self.__validate)
        self.entry.grid(sticky = tk.E + tk.W)
        self.scale = tk.Scale(self, from_ = value, to = value, command = self.__update,
                orient = tk.HORIZONTAL, resolution = 0.1)
        self.scale.grid(sticky = tk.E + tk.W)
        self.min = value
        self.max = value
        self.values = [value]
        self.get = self.scale.get
    def __update(self, value):
        if callable(self.command):
            self.command(self.name, value)
    def __validate(self, value):
        self.values = np.array(value.split(), dtype = float)
        self.min = min(self.values)
        self.max = max(self.values)
        self.scale.configure(from_ = self.min, to = self.max)

class LSelect(matplotlib.widgets.Lasso):
    def onrelease(self, event):
        if self.ignore(event):
            return
        if self.verts is not None:
            self.verts.append((event.xdata, event.ydata))
            self.callback(self.verts)
            self.ax.lines.remove(self.line)
        self.verts = None
        self.disconnect_events()

class DSFrame(tk.LabelFrame):
    def __init__(self, master, system, command = None, **kwargs):
        tk.LabelFrame.__init__(self, master, text = "System", **kwargs)
        self.columnconfigure(0, weight=1)
        self.system = system
        self.command = command
        self.eqn_x = tk.StringVar(self, self.system.get(1))
        self.eqn_y = tk.StringVar(self, self.system.get(2))
        self.params = dict()
        self.preset = tk.StringVar(self)
        choices = list(dynsystem.presets.keys())
        choices.sort()
        tk.OptionMenu(self, self.preset, *choices,
                command = self._load_preset).grid(sticky = tk.E+ tk.W, columnspan = 2)
        tk.Label(self, text = "x_dot:").grid()
        self.entry_x = VEntry(self, self.eqn_x, validator = self._update_system)
        self.entry_x.grid(row=1, column = 1)
        tk.Label(self, text = "y_dot:").grid()
        self.entry_y = VEntry(self, self.eqn_y, validator = self._update_system)
        self.entry_y.grid(row=2, column = 1)
        self._load_preset(choices[-1])

    def _update_system(self, *args):
        x_dot = self.eqn_x.get()
        y_dot = self.eqn_y.get()
        self.system.update(x_dot, y_dot)
        self._update_params()

    def _update_system_params(self, *args):
        for p in self.params:
            self.system.params[p] = self.params[p].get()
        if callable(self.command):
            self.command()

    def _load_preset(self, name):
        self.preset.set(name)
        x_dot, y_dot, params = dynsystem.presets[name]
        self.system.params.update(params)
        self.eqn_x.set(x_dot)
        self.eqn_y.set(y_dot)
        for p in params:
            self.params[p].entry.set(params[p])

    def _update_params(self, *args):
        for p in list(self.params.keys()):
            if not p in self.system.params:
                self.params[p].destroy()
                del self.params[p]
        for p in sorted(self.system.params):
            if not p in self.params:
                pe = Param(self, p, self.system.params[p], command = self._update_system_params)
                pe.grid(columnspan = 2, sticky = tk.E+tk.W)
                self.params[p] = pe

def test_system():
    x = e1.get()
    y = e2.get()
    u, v = d([x, y])
    print("x: %g, y: %g, u: %g, v:%g" % (x, y, u, v))
    print(d.params)

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
    dsf_test()
