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
        self.anim_timer   = self.fig.canvas.new_timer(interval = 5)
        self.anim_info    = tk.StringVar(self.root, "")
        self.anim_running = tk.BooleanVar(self.root, False)
        self.anim_tstep   = tk.IntVar(self.root, 0)
        self.anim_time    = tk.DoubleVar(self.root, 0.0)
        self.anim_tmax    = 0.0
        self.anim_timer.add_callback(self.anim_update)
        self.locations = []
        self.root.columnconfigure(0, weight=1)
        self.fig.bind('button_press_event', self.button_press)
        self.fig.bind('button_release_event', self.button_release)
        self.fig.bind('scroll_event', self.button_scroll)
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
        opts.polar          = tk.BooleanVar(self.root, False)
        opts.domain_factor  = tk.DoubleVar(self.root, 1.5)
        return opts

    def _computational_domain(self):
        graph_limits = self.fig.get_limits()
        comp_limits = []
        factor = self.opts.domain_factor.get()
        for l in graph_limits:
            comp_limits.append(helpers.scale_domain(l, factor))
        return comp_limits

    def update_trajectories(self, *args):
        picked = self.locations
        self.locations = []
        self.fig.trajectories = []
        self.anim_tmax = 0.0
        for pos in picked:
            self.handle_pick(pos)
        self.update_fig()

    def update_fig(self, *args):
        self.fig.clear()
        (x1, x2), (y1, y2) = self._computational_domain()[:2]
        x, y = np.meshgrid(np.linspace(x1, x2, 100), np.linspace(y1, y2, 100))
        u, v = self.system((x,y))
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
        self.fig.clear()
        self.update_fig()

    def _set_proj(self, *args):
        if self.opts.polar.get():
            self.fig.set_proj('polar')
        else:
            self.fig.set_proj('rectilinear')
        self._reset_fig()

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

    def handle_pick(self, pos):
        threshold = 1e-4
        self.locations.append(pos)
        traj = self.system.trajectory(pos, self.opts.tmax.get(), threshold = threshold,
                limits = self._computational_domain(), bidirectional = True)
        if traj.dist[-1] < 10*threshold:
            return
        if traj.t[-1] > self.anim_tmax:
            self.anim_tmax = traj.t[-1]
        self.fig.add_trajectory(traj, style = 'r-')
        self.fig.draw_trajectory(traj)
        self.fig.draw()

    def button_scroll(self, evt):
        pass
    def button_press(self, evt):
        pass
    def button_release(self, evt):
        # >HACK< Ignore if we are panning/zooming using the toolbar.
        if evt.button != 1:
            return
        #if self.fig.canvas.toolbar.mode != '':
        #    return
        if not evt.inaxes is self.fig.ax_main:
            return
        pos = evt.xdata, evt.ydata
        self.handle_pick(pos)

    def _add_widgets(self, frame):
        f_system = DSFrame(frame, self.system, command = self.update_trajectories)
        f_system.grid(sticky = tk.W + tk.E)

        f_controls = tk.Frame(frame)
        f_controls.grid(sticky = tk.W + tk.E)
        tk.Checkbutton(f_controls, text = "Nullclines", variable = self.opts.nullclines,
                command = self.update_fig).grid(row=0, column=0)
        tk.Checkbutton(f_controls, text = "Quiver", variable = self.opts.quiver,
                command = self.update_fig).grid(row = 0, column = 1)
        tk.Checkbutton(f_controls, text = "Polar", variable = self.opts.polar,
                command = self._set_proj).grid(columnspan = 2)
        tk.Button(f_controls, text = "Reset", underline = 0,
                command = self._reset_fig).grid(row=1,column=2)

        # Trajectories frame.
        row = 0
        f_traj = tk.LabelFrame(frame, text = "Trajectories:")
        f_traj.grid(sticky = tk.E+tk.W)
        f_traj.columnconfigure(2,weight=1)
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
