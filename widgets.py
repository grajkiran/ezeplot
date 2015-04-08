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
import presets

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
            self.errmsg.set("")
        except BaseException as e:
            self['bg'] = '#ff5555'
            self.errmsg.set(str(e))
            if self.debug:
                sys.stderr.write("%s\n" % self.errmsg.get())

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

class Param(tk.Frame):
    def __init__(self, master, name, value = 1.0, validator = None,
            command = None):
        tk.Frame.__init__(self, master)
        self.name = name
        self.var = tk.DoubleVar(self, value)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.label = tk.Label(self, text = name, width = 6)
        self.entry = VEntry(self, self.var, width = 6, validator = validator,
                command = command)
        self.label.grid(sticky = tk.E + tk.W)
        self.entry.grid(row = 0, column = 1, sticky = tk.E+tk.W)
        self.get = self.var.get
        self.set = self.var.set
    def enable(self, name = None, value = None):
        self.label.configure(state = tk.NORMAL)
        if name is None:
            name = self.name
        self.label.configure(text = name)
        self.entry.configure(state = tk.NORMAL)
        if value is not None:
            self.var.set(value)
    def disable(self):
        self.label.configure(state = tk.DISABLED, text = "")
        self.entry.configure(state = tk.DISABLED)

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
    def __init__(self, master, system, command = None, n_params = 10, preset_cmd = None, **kwargs):
        tk.LabelFrame.__init__(self, master, text = "System", **kwargs)
        self.columnconfigure(1, weight=1)
        self.system = system
        self.command = command
        self.preset_cmd = preset_cmd
        self.eqn_x = tk.StringVar(self, self.system.get(1))
        self.eqn_y = tk.StringVar(self, self.system.get(2))
        self.eqn_z = tk.StringVar(self, self.system.get(3))
        self.preset = tk.StringVar(self)
        choices = list(presets.systems.keys())
        choices.sort()
        choices.append("User defined")
        tk.OptionMenu(self, self.preset, *choices,
                command = self._load_preset).grid(sticky = tk.E+ tk.W, columnspan = 2)
        tk.Label(self, text = "x_dot:").grid()
        self.entry_x = VEntry(self, self.eqn_x, command = self.command,
                validator = self._update_system)
        self.entry_x.grid(row=1, column = 1, sticky = tk.E + tk.W)
        tk.Label(self, text = "y_dot:").grid()
        self.entry_y = VEntry(self, self.eqn_y, command = self.command,
                validator = self._update_system)
        self.entry_y.grid(row=2, column = 1, sticky = tk.E + tk.W)
        tk.Label(self, text = "z_dot:").grid()
        self.entry_z = VEntry(self, self.eqn_z, command = self.command,
                validator = self._update_system)
        self.entry_z.grid(row=3, column = 1, sticky = tk.E + tk.W)
        self.params = list()
        pframe = tk.Frame(self)
        pframe.grid(columnspan = 2, sticky = tk.E + tk.W)
        pframe.columnconfigure(0, weight=1)
        pframe.columnconfigure(1, weight=1)
        for row in range(n_params//2):
            for col in 0, 1:
                pe = Param(pframe, "", 0.0, command = self.command,
                        validator = self._update_system_params)
                pe.grid(row = row, column = col, sticky = tk.E + tk.W)
                self.params.append(pe)
        self._load_preset(choices[-1])

    def _load_preset(self, name):
        self.preset.set(name)
        if name != "User defined":
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
            name = pe.label.cget('text')
            value = pe.get()
            if name in self.system.params:
                self.system.params[name] = value

    def _update_params(self, *args):
        items = sorted(self.system.params.items())
        for i in range(len(items)):
            self.params[i].enable(name = items[i][0], value = items[i][1])
        for pe in self.params[len(items):]:
            pe.disable()

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
