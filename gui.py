#!/usr/bin/python3
import sys

try:
    import tkinter as tk
except:
    import Tkinter as tk
import numpy as np
from widgets import *

class Options(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__

class AppWindow():
    def __init__(self, root, system, fig):
        self.system = system
        self.root = root
        self.fig = fig
        self.opts = self._init_options()
        self.anim_timer   = self.fig.canvas.new_timer(interval = 5)
        self.anim_info    = tk.StringVar(self.root, "")
        self.anim_running = tk.BooleanVar(self.root, False)
        self.anim_tstep   = tk.IntVar(self.root, 0)
        self.anim_time    = tk.DoubleVar(self.root, 0.0)
        self.anim_timer.add_callback(self.anim_update)
        self.root.columnconfigure(0, weight=1)
        self.fig.bind('button_release_event', self.handle_mouse)
        self._add_widgets(self.root)

    def _init_options(self, fname = None):
        opts = Options()
        opts.tmax           = tk.DoubleVar(self.root, 5)
        opts.dt             = tk.DoubleVar(self.root, 0.05)
        opts.coarsen        = tk.IntVar(self.root, 1)
        opts.quiver         = tk.BooleanVar(self.root, False)
        opts.quiver_scale   = tk.BooleanVar(self.root, True)
        opts.nullclines     = tk.BooleanVar(self.root, True)
        opts.spacing        = tk.DoubleVar(self.root, 0.3)
        opts.temporal       = tk.BooleanVar(self.root, False)
        opts.polar          = tk.BooleanVar(self.root, False)
        return opts

    def nullclines(self, *args):
        field = None
        if self.opts.nullclines.get():
            field = compute_field()
        fig.nullclines(field)

    def quiver(self, *args):
        field = None
        if self.opts.quiver.get():
            field = compute_field()
        fig.quiver(field, self.opts.quiver_scale.get())

    def _set_mode(*args): pass
    def _reset_plot(*args): pass
    def update_trajectories(*args): pass

    def anim_update(self):
        tstep = self.anim_tstep.get()
        anim_time = tstep * self.opts.dt.get()
        if anim_time > self.opts.tmax.get():
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

    def handle_mouse(self, evt):
        # >HACK< Ignore if we are panning/zooming using the toolbar.
        if evt.button != 1:
            return
        if self.fig.canvas.toolbar.mode != '':
            return
        if not evt.inaxes is self.fig.ax_main:
            return
        pos = evt.xdata, evt.ydata
        traj1 = self.system.trajectory(pos, 50, threshold = 1e-4,
                limits = self.fig.get_limits(), bidirectional = True)
        self.fig.add_trajectory(traj1, style = 'r-')

    def _add_widgets(self, frame):
        f_system = DSFrame(frame, self.system, command = self.update_trajectories)
        f_system.grid(sticky = tk.W + tk.E)

        f_controls = tk.Frame(frame)
        f_controls.grid(sticky = tk.W + tk.E)
        tk.Checkbutton(f_controls, text = "Nullclines", variable = self.opts.nullclines,
                command = self.nullclines).grid(row=0, column=0)
        tk.Checkbutton(f_controls, text = "Quiver", variable = self.opts.quiver,
                command = self.quiver).grid(row = 0, column = 1)
        tk.Checkbutton(f_controls, text = "Uniform",
                variable = self.opts.quiver_scale).grid(row = 0, column = 2)
        tk.Checkbutton(f_controls, text = "Polar mode", variable = self.opts.polar,
                command = self._set_mode).grid(columnspan = 2)
        tk.Button(f_controls, text = "Reset", underline = 0,
                command = self._reset_plot).grid(row=1,column=2)

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

class VectorFieldOLD:
    def __init__(self, system, xlim = (-1.0, 1.0), ylim = (-1.0, 1.0)):
        self.system = system
        self.trajectories = []

        self.fig = plt.figure()
        self.fig.canvas.set_window_title("2D Autonomous Systems")
        self._blit = self.fig.canvas.supports_blit
        tk_master = self.fig.canvas.get_tk_widget()
        self.anim_timer = self.fig.canvas.new_timer(interval = 5)
        self.anim_timer.add_callback(self.anim_update)
        self.anim_curr_step = 0
        self.anim_dt = tk.DoubleVar(tk_master, value = 0.05)
        self.anim_tsteps = tk.IntVar(tk_master, value = 500)
        self.anim_cmd = tk.StringVar(tk_master, value = "Animate")
        self.anim_info = tk.StringVar(tk_master, value = "")
        self.anim_running = False
        self.anim_coarsen = 5
        self.quiver = tk.IntVar(tk_master, value = 0)
        self.quiver_scale = tk.IntVar(tk_master, value = 1)
        self.nullclines = tk.IntVar(tk_master, value = 0)
        self.auto_plot = tk.IntVar(tk_master, value = 1)
        self.spacing = tk.DoubleVar(tk_master, value = 0.3)
        self.temporal = tk.IntVar(tk_master, value = 0)
        self.polar = tk.IntVar(tk_master, value = 0)

        self.ax_vis = self.fig.add_subplot(222, label = "Visual", visible = False)
        self.ax_tx = self.fig.add_subplot(223, label = "TX plane",
                title = "TX plane", visible = False)
        self.ax_ty = self.fig.add_subplot(224, label = "TY plane",
                title = "TY plane", visible = False, sharex = self.ax_tx)
        self.ax_rect = self.fig.add_subplot(221, xlim = xlim, ylim = ylim,
                label = "Phase Plane", title = "Phase Plane")
        self.ax_polar = self.fig.add_subplot(221, xlim = xlim, ylim = ylim, polar = True,
                label = "Phase Plane", title = "Phase Plane")
        self.ax = self.ax_rect
        self.fig.subplots_adjust(0.06, 0.01, 0.99, 0.94, 0.1, 0.1)
        self._set_mode()
        self._clear_axes(temporal = True)
        self._add_widgets()
        self.fig.canvas.mpl_connect("button_press_event", self.handle_mouse)
        self.fig.canvas.mpl_connect("scroll_event", self.handle_scroll)

    def handle_scroll(self, evt):
        def extend_limits(lim, ratio = 1.0):
            a, b = lim
            center = (a+b)/2.0
            d1 = ratio * (center - a)
            d2 = ratio * (b - center)
            return center - d1, d2 + center
        if evt.button == "up":
            scale = 0.8
        else:
            scale = 1.25
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        if self.polar.get():
            self.ax.set_ylim(0, scale*ylim[1])
        else:
            self.ax.set_xlim(extend_limits(xlim, scale))
            self.ax.set_ylim(extend_limits(ylim, scale))
        self.fig.canvas.draw()

    def update_trajectories(self):
        for i, t in enumerate(self.trajectories):
            traj = self.initialize_trajectory(t.start, direction = t.direction)
            self.trajectories[i] = traj
        self.plot_trajectories()

    def anim_update(self):
        tstep = self.anim_curr_step
        self.anim_curr_step += 1
        if self.anim_curr_step > self.anim_tsteps.get():
            self.anim_curr_step -= self.anim_tsteps.get()
        self.anim_info.set("%d" % tstep)
        if self._blit:
            self.fig.canvas.restore_region(self.blank_bg)
        else:
            if tstep % self.anim_coarsen != 0: return
            self._clear_axes(temporal = True, draw = False)
        for traj in self.trajectories:
            self.draw_trajectory(traj, tstep, arrows = 0, marker = True)
            if self.temporal.get():
                if self._blit:
                    self.fig.canvas.restore_region(self.blank_tx)
                    self.fig.canvas.restore_region(self.blank_ty)
                self.draw_temporal(traj, tstep, marker = True)
        self.blit()

    def toggle_traj_animation(self):
        if self.anim_cmd.get() == 'Animate':
            self.anim_timer.start()
            self.anim_cmd.set('Stop')
        else:
            self.anim_timer.stop()
            self.anim_cmd.set('Animate')

#    def toggle_traj_animation_old(self):
#        if self.anim_running:
#            self.stop_traj_animation()
#        else:
#            self.anim_thread = threading.Thread(target = self.animated_trajectories)
#            self.stop_animation = False
#            self.anim_thread.start()
#    def stop_traj_animation(self):
#        self.stop_animation = True
#        time.sleep(0.05)
#
#    def animated_trajectories(self):
#        self._clear_axes()
#        self.anim_running = True
#        self.anim_cmd.set("Stop")
#        self.stop_animation = False
#        remaining = len(self.trajectories)
#        for tstep in range(1, 2*self.anim_tsteps.get()):
#            self.anim_info.set("%d" % tstep)
#            if remaining <= 0 or self.stop_animation:
#                break
#            if self._blit:
#                self.fig.canvas.restore_region(self.blank_bg)
#            else:
#                if tstep % self.anim_coarsen != 0: continue
#                self._clear_axes(temporal = True, draw = False)
#            for traj in self.trajectories:
#                if len(traj.x) == tstep:
#                    remaining -= 1
#                self.draw_trajectory(traj, tstep, arrows = 0, marker = True)
#                if self.temporal.get():
#                    if self._blit:
#                        self.fig.canvas.restore_region(self.blank_tx)
#                        self.fig.canvas.restore_region(self.blank_ty)
#                    self.draw_temporal(traj, tstep, marker = True)
#            self.blit()
#        self.anim_info.set("")
#        self.anim_running = False
#        self.anim_cmd.set("Animate")

    def initialize_trajectory(self, loc, direction = 1, **kwargs):
        """Compute the trajectory, initialize the trail and head."""
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        if self.polar.get():
            xlim = None
        traj = self.system.trajectory(loc, dt = self.anim_dt.get(),
                xlim = xlim, ylim = ylim,
                max_steps = self.anim_tsteps.get(), direction = direction)
        return traj

    def draw_artist(self, artist, ax = None):
        if self._blit:
            if ax is None:
                self.ax.draw_artist(artist)
            else:
                ax.draw_artist(artist)

    def blit(self):
        if self._blit:
            self.fig.canvas.blit()
        else:
            #self.fig.canvas.draw()
            self.fig.canvas.draw_idle()

    def draw_temporal(self, traj, n = 0, marker = False):
        if not self.temporal.get():
            return
        if n == 0 or n >= len(traj.x):
            n = len(traj.x)
        lx, = self.ax_tx.plot(traj.t[:n], traj.x[:n], 'r-')
        ly, = self.ax_ty.plot(traj.t[:n], traj.y[:n], 'r-')
        if self._blit:
            self.ax_tx.draw_artist(lx)
            self.ax_ty.draw_artist(ly)
        if marker:
            mx, = self.ax_tx.plot(traj.t[n-1:n], traj.x[n-1:n], 'ro')
            my, = self.ax_ty.plot(traj.t[n-1:n], traj.y[n-1:n], 'ro')
            if self._blit:
                self.ax_tx.draw_artist(mx)
                self.ax_ty.draw_artist(my)

    def draw_trajectory(self, traj, n = 0, marker = False, arrows = 0, trail = True):
        if n == 0 or n >= len(traj.x):
            n = len(traj.x)
        x = traj.x[:n]
        y = traj.y[:n]
        if trail:
            traj.trail, = self.ax.plot(x, y, 'r-')
            self.draw_artist(traj.trail)
        if marker:
            traj.marker, = self.ax.plot(x[-1:], y[-1:], 'ro', ms = 7)
            self.draw_artist(traj.marker)
        if arrows == 0:
            return
        # Compute arrow positions
        if arrows > 1:
            if arrows >= n:
                arrows = n-1
            tails = curve_resample(x, y, count = arrows+2)[0][1:-1]
            heads = tails + 1
        else:
            tails = [n-1]
            heads = [n]
        # Draw arrows
        for i in tails:
            from_ = x[i], y[i]
            to_ = traj.x[i+1], traj.y[i+1]
            dx = 0.01 * (to_[0] - from_[0])
            dy = 0.01 * (to_[1] - from_[1])
            to_ = from_[0] + dx, from_[1] + dy
            artist = self.ax.annotate("", xy = to_, xytext = from_,
                    arrowprops = dict(arrowstyle = "->"))
            self.draw_artist(artist)

    def compute_field(self, n1 = None, n2 = None, density = 1.0):
        x1, x2 = self.ax.get_xlim()
        y1, y2 = self.ax.get_ylim()
        n1 = int(density * (n1 or abs(x1-x2)/self.spacing.get()))
        n2 = int(density * (n2 or abs(y1-y2)/self.spacing.get()))
        mesh = np.meshgrid(np.linspace(x1, x2, n1), np.linspace(y1, y2, n2))
        x, y = mesh
        u, v = self.system(mesh)
        return x, y, u, v

    def add_nullclines(self, *args):
        if self.nullclines.get() == 0:
            return
        x, y, u, v = self.compute_field(density = 5)
        nc_u = self.ax.contour(x, y, u, [0.0], linewidths = 2, colors = 'b')#.collections[0].get_paths()[0].vertices
        nc_v = self.ax.contour(x, y, v, [0.0], linewidths = 2, colors = 'r')#.collections[0].get_paths()[0].vertices
        self.draw_artist(nc_u.collections[0])
        self.draw_artist(nc_v.collections[0])
        self.blit()

    def add_quiver(self, mesh = None):
        if self.quiver.get() == 0:
            return
        if mesh is None:
            x, y, u, v = self.compute_field(25, 25)
        else:
            x, y = mesh
            u, v = self.system(mesh)
        if self.quiver_scale.get():
            mag = np.ma.masked_less(np.sqrt(u*u+v*v), 1.0e-4)
            u = u/mag
            v = v/mag
        artist = self.ax.quiver(x,y,u,v, units = 'xy')
        self.draw_artist(artist)
        self.blit()

    def random(self, count):
        count = int(count)
        x1, x2 = self.ax.get_xlim()
        y1, y2 = self.ax.get_ylim()
        x = np.random.random(count) * (x2 - x1) + x1
        y = np.random.random(count) * (y2 - y1) + y1
        return np.transpose(np.array([x, y]))

    def event_info(self, event = None):
        print(self.anim_toggle.get(), self.anim_status.get())
        print(event)

def test_vectorfield():
    x1, x2 = (-5.0, 5.0)
    y1, y2 = (-5.0, 5.0)
    d = dynsystem.DynamicSystem(x_dot = "x * (3 - x - y)", y_dot = "y * (x - 1)")
    f = VectorField(d, xlim = (x1, x2), ylim = (y1, y2))
    f.fig.show()
    f.root.mainloop()
    return d, f

if __name__ == '__main__':
    np.seterr(invalid = 'print')
    d, f = test_vectorfield()
