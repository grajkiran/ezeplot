#!/usr/bin/python3
import sys

try:
    import tkinter as tk
except:
    import Tkinter as tk
import numpy as np
from widgets import *
import helpers
import plotting

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
        #self.limits_dialog = PlotLimits(self.root, self.fig, self.opts.limits)
        self.anim_timer   = self.fig.canvas.new_timer(interval = 5)
        self.anim_info    = tk.StringVar(self.root, "")
        self.anim_running = tk.BooleanVar(self.root, False)
        self.anim_tstep   = tk.IntVar(self.root, 0)
        self.anim_time    = tk.DoubleVar(self.root, 0.0)
        self.anim_tmax    = 0.0
        self.anim_timer.add_callback(self.anim_update)
        self.locations = []
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
        self.update_fig()

    def _init_options(self, fname = None):
        opts = Options()
        opts.tmax           = tk.DoubleVar(self.root, 25)
        opts.dt             = tk.DoubleVar(self.root, 0.05)
        opts.coarsen        = tk.IntVar(self.root, 1)
        opts.quiver         = tk.BooleanVar(self.root, False)
        opts.quiver_scale   = tk.BooleanVar(self.root, True)
        opts.nullclines     = tk.BooleanVar(self.root, True)
        opts.spacing        = tk.DoubleVar(self.root, 0.3)
        opts.temporal       = tk.BooleanVar(self.root, False)
        opts.projection     = tk.StringVar(self.root, '2D')
        #opts.domain_factor  = tk.DoubleVar(self.root, 1.5)
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

    def _update_system_limits(self, *args, prompt = True):
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

#    def _computational_domain(self):
#        graph_limits = self.fig.get_limits()
#        comp_limits = []
#        factor = self.opts.limits.factor.get()
#        for l in graph_limits:
#            comp_limits.append(helpers.scale_domain(l, factor))
#        self.system.limits = comp_limits
#        return comp_limits

    def update_trajectories(self, *args):
        picked = self.locations
        self.locations = []
        self.fig.trajectories = []
        self.anim_tmax = 0.0
        for pos in picked:
            self.add_location(pos)
        self.update_fig()

    def update_fig(self, *args):
        self.fig.clear(tmax = self.opts.tmax.get())
        #(x1, x2), (y1, y2) = self._computational_domain()[:2]
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
        for t in self.fig.trajectories:
            self.fig.draw_trajectory(t)
        self.fig.draw()

    def _reset_fig(self, *args):
        self.anim_running.set(False)
        self.anim_timer.stop()
        self.locations = []
        self.fig.trajectories = []
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
        if anim_time > self.anim_tmax:
            anim_time = 0.0
            self.anim_tstep.set(0)
        else:
            self.anim_tstep.set(tstep+1)
        self.anim_info.set("%0.3f" % anim_time)
        self.fig.anim_update(anim_time)

    def toggle_traj_animation(self):
        if not self.anim_running.get():
            self.anim_timer.stop()
        else:
            self.anim_timer.start()

    def add_location(self, pos = None):
        threshold = 1e-4
        if pos is None or isinstance(pos, tk.Event):
            pos = helpers.parse_coords(self.location_str.get())
        if len(pos) == 2:
            pos = pos[0], pos[1], 0.0
        if pos in self.locations:
            print("Trajectory already exists.")
            return
#        limits = self._computational_domain()
#        if self.opts.projection.get() == 'Polar':
#            limits[0] = None
        try:
            traj = self.system.trajectory(pos, self.opts.tmax.get(), threshold = threshold,
                bidirectional = True, nsteps = 5 * (self.opts.tmax.get()/self.opts.dt.get()),
                max_step = self.opts.dt.get())
        except:
            pos_str = ", ".join(map(str, pos))
            sys.stderr.write("Could not compute trajectory from: %s\n" % pos_str)
            print(sys.exc_info())
            return
        self.locations.append(pos)
        if traj.dist[-1] < 10*threshold:
            return
        if traj.t[-1] > self.anim_tmax:
            self.anim_tmax = traj.t[-1]
        self.fig.add_trajectory(traj, style = 'r-')
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
        return
        #print(evt.inaxes, evt.x, evt.y, evt.xdata, evt.ydata, evt.button)
    def button_release(self, evt):
        if evt.button != 1:
            return
        if not evt.inaxes is self.fig.ax_main:
            return
        if self.fig.ax_main is self.fig.ax_3d:
            return
        if self.mouse_mode == 'pick':
            pos = evt.xdata, evt.ydata
            self.add_location(pos)
        else:
            self.mouse_mode = 'pick'

    def _add_widgets(self, frame):
        f_system = DSFrame(frame, self.system, command = self.update_trajectories)
        f_system.grid(sticky = tk.W + tk.E)

        f_controls = tk.Frame(frame)
        f_controls.grid(sticky = tk.W + tk.E)
        tk.Checkbutton(f_controls, text = "Nullclines", variable = self.opts.nullclines,
                command = self.update_fig).grid(row=0, column=0)
        tk.Checkbutton(f_controls, text = "Quiver", variable = self.opts.quiver,
                command = self.update_fig).grid(row = 0, column = 1)
        tk.Checkbutton(f_controls, text = "Temporal", variable = self.opts.temporal,
                command = self._set_temporal).grid(row = 0, column = 2)
        optmenu = tk.OptionMenu(f_controls, self.opts.projection,
                *PROJECTIONS.keys(), command = self._set_proj)
        optmenu.configure(width = 5)
        optmenu.grid(columnspan = 2)
        tk.Button(f_controls, text = "Reset", underline = 0,
                command = self._reset_fig).grid(row=1,column=2)
        tk.Button(f_controls, text = "Change plot limits",
                command = self._update_system_limits).grid(columnspan = 3)

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
        tk.Label(f_traj, textvariable = self.anim_info).grid(row = row, column = 2)
        row += 1
        tk.Button(f_traj, text = "Update", command = self.update_trajectories).grid(row = row, column = 0)
        tk.Button(f_traj, text = "Restart", command = lambda: self.anim_tstep.set(0)).grid(row = row, column = 1)
        tk.Checkbutton(f_traj, text = "Animated", variable = self.anim_running,
                command = self.toggle_traj_animation).grid(sticky = tk.E + tk.W, row = row, column = 2)

if __name__ == '__main__':
    np.seterr(invalid = 'print')
