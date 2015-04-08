#!/usr/bin/python3
import sys

try:
    import tkinter as tk
    from tkinter import ttk
except:
    import Tkinter as tk
    import ttk
import numpy as np
from widgets import *
import helpers
import plotting
from poincare import PWindow
import presets

PROJECTIONS = dict({'2D': 'rect', 'Polar': 'polar', '3D': '3d'})

class Options(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__

class AppWindow():
    def __init__(self, root, system, blit = True, embedded = False):
        self.system = system
        self.root = root
        master = self.root if embedded else None
        self.fig = plotting.Figure(blit = blit, master = master)
        self.opts = self._init_options()
        self._update_system_limits(prompt = False)
        self.anim_timer   = self.fig.canvas.new_timer(interval = 5)
        self.anim_info    = tk.StringVar(self.root, "")
        self.anim_running = tk.BooleanVar(self.root, False)
        self.anim_tstep   = tk.IntVar(self.root, 0)
        self.anim_tmax    = 0.0
        self.anim_timer.add_callback(self.anim_update)
        self.pointer_info = tk.StringVar(self.root, "")
        self.trajectories = dict()
        # Used for manually adding a point.
        self.location_str = tk.StringVar(self.root, "")
        self.mouse_mode = "pick" # pan, dynamic
        self.root.columnconfigure(0, weight=1)
        self.fig.bind('button_press_event', self.button_press)
        self.fig.bind('button_release_event', self.button_release)
        self.fig.bind('motion_notify_event', self.mouse_move)
        self.fig.bind('resize_event', self.update_fig)
        if embedded:
            cframe = tk.Frame(self.root)
            cframe.pack()
            self._add_widgets(cframe)
        else:
            self._add_widgets(self.root)
            self.fig.bind('close_event', lambda evt: root.quit())
        tk.Label(self.root, textvariable = self.pointer_info, anchor = tk.W,
                relief = tk.SUNKEN,).pack(side = tk.BOTTOM, fill = tk.X)
        self.update_fig()

    def _init_options(self, fname = None):
        opts = Options()
        opts.tmax           = tk.DoubleVar(self.root, 25)
        opts.dt             = tk.DoubleVar(self.root, 0.05)
        opts.quiver         = tk.BooleanVar(self.root, False)
        opts.nullclines     = tk.BooleanVar(self.root, False)
        opts.temporal       = tk.BooleanVar(self.root, False)
        opts.projection     = tk.StringVar(self.root, '2D')
        opts.limits         = Options()
        opts.limits.factor  = tk.DoubleVar(self.root, 1.5)
        opts.limits.xmin    = tk.DoubleVar(self.root, -5.0)
        opts.limits.xmax    = tk.DoubleVar(self.root, 5.0)
        opts.limits.ymin    = tk.DoubleVar(self.root, -5.0)
        opts.limits.ymax    = tk.DoubleVar(self.root, 5.0)
        opts.limits.zmin    = tk.DoubleVar(self.root, -5.0)
        opts.limits.zmax    = tk.DoubleVar(self.root, 5.0)
        opts.limits.per_x   = tk.BooleanVar(self.root, False)
        opts.limits.per_y   = tk.BooleanVar(self.root, False)
        opts.limits.per_z   = tk.BooleanVar(self.root, False)
        return opts

    def _load_preset(self, name):
        opts = self.opts
        preset = presets.systems[name]
        for opt in ('tmax', 'dt', 'projection'):
            if opt in preset:
                opts[opt].set(preset[opt])
        if 'xlim' in preset:
            x1, x2 = preset['xlim']
            opts.limits.xmin.set(x1)
            opts.limits.xmax.set(x2)
        if 'ylim' in preset:
            y1, y2 = preset['ylim']
            opts.limits.ymin.set(y1)
            opts.limits.ymax.set(y2)
        if 'zlim' in preset:
            z1, z2 = preset['zlim']
            opts.limits.zmin.set(z1)
            opts.limits.zmax.set(z2)
        #self._update_system_limits(prompt = False)
        self._reset_fig()
        self._set_proj()
        if 'locations' in preset:
            for pos in preset['locations']:
                self.add_location(pos)
        #self.update_trajectories()

    def _update_system_limits(self, evt = None, prompt = True):
        if prompt:
            limits_dialog = PlotLimits(self.root, self.fig, self.opts.limits)
        limits = self.opts.limits
        xmin = limits.xmin.get()
        xmax = limits.xmax.get()
        ymin = limits.ymin.get()
        ymax = limits.ymax.get()
        zmin = limits.zmin.get()
        zmax = limits.zmax.get()
        factor = limits.factor.get()
        plot_limits = [(xmin, xmax), (ymin, ymax), (zmin, zmax)]
        periodic = [limits.per_x.get(), limits.per_y.get(), limits.per_z.get()]
        self.fig.set_limits(*plot_limits)
        self.fig.draw(force = True)
        for i in range(3):
            if periodic[i]:
                self.system.limits[i] = plot_limits[i]
            else:
                self.system.limits[i] = helpers.scale_domain(plot_limits[i], factor)
        self.system.periodic = periodic
        if self.opts.projection.get() == 'Polar':
            # Ignore xlim (theta) in polar mode.
            self.system.limits[0] = None
            self.system.periodic[0] = True

    def show_poincare_dialog(self, *args):
        l = self.opts.limits
        x1, x2 = l.xmin.get(), l.xmax.get()
        y1, y2 = l.ymin.get(), l.ymax.get()
        z1, z2 = l.zmin.get(), l.zmax.get()
        w = PWindow(self.root, self.trajectories,
                ((x1, x2), (y1, y2), (z1, z2)))

    def update_trajectories(self, *args):
        picked = list(self.trajectories.keys())
        self.trajectories.clear()
        self.anim_tmax = 0.0
        self.update_fig()
        for pos in picked:
            self.add_location(pos)

    def update_fig(self, *args):
        self.fig.clear(tmax = self.opts.tmax.get())
        x1 = self.opts.limits.xmin.get()
        x2 = self.opts.limits.xmax.get()
        y1 = self.opts.limits.ymin.get()
        y2 = self.opts.limits.ymax.get()
        x, y = np.meshgrid(np.linspace(x1, x2, 100), np.linspace(y1, y2, 100))
        u, v, w = self.system((x,y))
        field = x, y, u, v
        field_quiv = [f[::4, ::4] for f in field]
        if self.opts.nullclines.get():
            self.fig.draw_nullclines(field)
        if self.opts.quiver.get():
            self.fig.draw_quiver(field_quiv)
        #NOTE: Add other bg elements like FP, LC, here before draw and save.
        self.fig.draw(force = True)
        self.fig.save()
        for t in self.trajectories.values():
            self.fig.draw_trajectory(t)
        self.fig.draw()

    def _reset_fig(self, *args):
        self.anim_running.set(False)
        self.anim_timer.stop()
        self.trajectories.clear()
        self.anim_tmax = 0.0
        self.fig.clear(tmax = self.opts.tmax.get())
        self.update_fig()

    def _set_temporal(self, *args):
        if self.opts.temporal.get():
            self.fig.set_mode(4)
        else:
            self.fig.set_mode(1)
        self.update_fig()

    def _set_proj(self, *args):
        proj = self.opts.projection.get()
        self.fig.set_proj(PROJECTIONS[proj])
        self.update_fig()
        self._update_system_limits(prompt = False)

    def anim_update(self):
        tstep = self.anim_tstep.get()
        anim_time = tstep * self.opts.dt.get()
        if anim_time > self.opts.tmax.get():
            anim_time = 0.0
            self.anim_tstep.set(0)
        else:
            self.anim_tstep.set(tstep+1)
        self.anim_info.set("%0.3f" % anim_time)
        self.fig.anim_update(anim_time, self.trajectories.values())

    def toggle_traj_animation(self):
        if not self.anim_running.get():
            self.anim_timer.stop()
            self.anim_tstep.set(0)
        else:
            self.anim_timer.start()

    def add_location(self, pos = None):
        styles = ["b", "g", "k", "m"]
        threshold = 1e-4
        if pos is None or isinstance(pos, tk.Event):
            pos = helpers.parse_coords(self.location_str.get())
        if len(pos) == 2:
            pos = pos[0], pos[1], 0.0
        pos = tuple(pos)
        if pos in self.trajectories:
            print("Trajectory already exists.")
            return
        try:
            traj = self.system.trajectory(pos, self.opts.tmax.get(), threshold = threshold,
                bidirectional = True, nsteps = 5 * (self.opts.tmax.get()/self.opts.dt.get()),
                max_step = self.opts.dt.get())
            #print("Computing trajectory took %g seconds" % t.seconds())
        except:
            pos_str = ", ".join(map(str, pos))
            sys.stderr.write("Could not compute trajectory from: %s\n" % pos_str)
            print(sys.exc_info())
            return
        if traj.dist[-1] < 10*threshold:
            return
        self.trajectories[pos] = traj
        if traj.t[-1] > self.anim_tmax:
            self.anim_tmax = traj.t[-1]
        style = (len(self.trajectories) - 1) % len(styles)
        self.fig.add_trajectory(traj, style = styles[style])
        self.fig.draw_trajectory(traj)
        self.fig.draw()

    def button_press(self, evt):
        if evt.key == 'shift' and evt.inaxes is self.fig.ax_main:
            self.mouse_mode = 'pan'
            self.pan_loc = evt.xdata, evt.ydata
        elif evt.key == 'control':
            self.mouse_mode = 'dynamic'
        else:
            self.mouse_mode = 'pick'
    def mouse_move(self, evt):
        s = ""
        if evt.inaxes is None:
            s = ""
        elif evt.inaxes == self.fig.ax_3d:
            s = ""
        elif evt.inaxes == self.fig.ax_main:
            s = evt.inaxes.format_coord(evt.xdata, evt.ydata)
        self.pointer_info.set(s)

    def button_release(self, evt):
        if evt.button != 1:
            return
        if not evt.inaxes is self.fig.ax_main:
            return
        if self.fig.ax_main is self.fig.ax_3d:
            self.update_fig()
            return
        if self.mouse_mode == 'pick':
            pos = evt.xdata, evt.ydata
            self.add_location(pos)
        else:
            self.mouse_mode = 'pick'

    def _add_widgets(self, frame):
        f_system = DSFrame(frame, self.system, preset_cmd = self._load_preset,
                command = self.update_trajectories)
        f_system.grid(sticky = tk.W + tk.E)

        f_controls = tk.Frame(frame)
        f_controls.grid(sticky = tk.W + tk.E)
        tk.Checkbutton(f_controls, text = "Nullclines", variable = self.opts.nullclines,
                command = self.update_fig).grid(row=0, column=0)
        tk.Checkbutton(f_controls, text = "Quiver", variable = self.opts.quiver,
                command = self.update_fig).grid(row = 0, column = 1)
        tk.Checkbutton(f_controls, text = "Graphs", variable = self.opts.temporal,
                command = self._set_temporal).grid(row = 0, column = 2)
        optmenu = tk.OptionMenu(f_controls, self.opts.projection,
                *PROJECTIONS.keys(), command = self._set_proj)
        optmenu.configure(width = 5)
        optmenu.grid(columnspan = 2)
        tk.Button(f_controls, text = "Limits",
                command = self._update_system_limits).grid()
        tk.Button(f_controls, text = "Poincare section",
                command = self.show_poincare_dialog).grid(row = 2, column = 1, columnspan = 2)

        # Trajectories frame.
        row = 0
        f_traj = tk.LabelFrame(frame, text = "Trajectories:")
        f_traj.grid(sticky = tk.E+tk.W)
        f_traj.columnconfigure(2,weight=1)
        tk.Label(f_traj, text = "Coords:").grid(row = row, column = 0)
        VEntry(f_traj, textvariable = self.location_str, validator = helpers.parse_coords,
                command = self.add_location, width = 12).grid(row = row, column = 1, columnspan = 2, sticky = tk.W)
        tk.Button(f_traj, text = "Add", command = self.add_location).grid(row = row, column = 2, sticky = tk.E)
        row += 1
        tk.Label(f_traj, text = "Time step:").grid(row = row, column = 0)
        VEntry(f_traj, textvariable = self.opts.dt, width = 5).grid(row = row, column = 1)
        row += 1
        tk.Label(f_traj, text = "Tmax:").grid(row = row, column = 0)
        VEntry(f_traj, textvariable = self.opts.tmax, width = 5).grid(row = row, column = 1)
        #tk.Label(f_traj, textvariable = self.anim_info).grid(row = row, column = 2)
        row += 1
        #tk.Button(f_traj, text = "Restart", command = lambda: self.anim_tstep.set(0)).grid(row = row, column = 1)
        tk.Checkbutton(f_traj, text = "Animate", pady = 2, padx = 4, font = "sans 12 bold", variable = self.anim_running, indicatoron = 0,
                command = self.toggle_traj_animation).grid(sticky = tk.E + tk.W, row = row, column = 2)

        f_update = tk.Frame(frame)
        f_update.columnconfigure(0, weight=1)
        tk.Button(f_update, text = "Update", command = self.update_trajectories,
                background = "#0000aa", activebackground = "#3333ff",
                foreground = "white", activeforeground = "white", font = "sans 16 bold",
                height = 1, width = 6).grid(
                        row = 0, columnspan=2, sticky = tk.S)
        tk.Button(f_update, text = "Reset", command = self._reset_fig, font = "sans 10 bold",
                background = "#aa0000", activebackground = "#ff5555",
                foreground = "white", activeforeground = "white").grid(
                        row = 0, column = 2, sticky = tk.E + tk.S)
        f_update.grid(sticky = tk.E + tk.W + tk.S)

if __name__ == '__main__':
    np.seterr(invalid = 'print')
