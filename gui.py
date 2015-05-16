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
import sys

try:
    import tkinter as tk
    from tkinter import ttk
    from tkinter.filedialog import asksaveasfilename
except:
    import Tkinter as tk
    import ttk
    from tkFileDialog import asksaveasfilename
import numpy as np
from widgets import *
import helpers
import plotting
from poincare import PWindow
import presets
import uptime
from timer import Timer
PROJECTIONS = dict({'2D': 'rect', u'Polar (x≡θ, y≡r)': 'polar', '3D': '3d'})

class Options(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__

class AppWindow():
    def __init__(self, root, system, blit = True):#, embedded = False):
        self.system = system
        self.root = root
        print(uptime.uptime(), "Creating figure window...")
        self.fig = plotting.Figure(root, blit = blit)
        print(uptime.uptime(), "Initializing options...")
        self.opts = self._init_options()
        self.anim_timer   = self.fig.canvas.new_timer(interval = 5)
        self.anim_info    = tk.StringVar(self.root, "")
        self.anim_running = tk.BooleanVar(self.root, False)
        self.anim_tstep   = tk.IntVar(self.root, 0)
        self.anim_tmax    = 0.0
        self.anim_timer.add_callback(self.anim_update)
        self.pointer_info = tk.StringVar(self.root, "")
        self.trajectories = dict()
        self.fixed_points = set()
        self.last_loc = 0.0, 0.0, 0.0
        self.press_loc = 0.0, 0.0
        # Used for manually adding a point.
        self.location_str = tk.StringVar(self.root, "")
        self.root.columnconfigure(0, weight=1)
        self.fig.bind('button_press_event', self.button_press)
        self.fig.bind('button_release_event', self.button_release)
        self.fig.bind('motion_notify_event', self.mouse_move)
        self.fig.bind('resize_event', self._handle_resize)
        pinfo_label = tk.Label(self.root, textvariable = self.pointer_info, anchor = tk.W,
                relief = tk.SUNKEN,)
        cframe = tk.Frame(self.root, bd = 0, relief = tk.RIDGE)
        self.menu = self._add_menubar()
        print(uptime.uptime(), "Creating widgets...")
        self.controls = self._add_widgets(cframe)
        print(uptime.uptime(), "Updating figure...")
        self._update_system_limits()
        self.update_fig()

        print(uptime.uptime(), "Finishing...")
        self.root.rowconfigure(1, weight = 1)
        self.root.columnconfigure(0, weight = 1)
        self.root.config(menu = self.menu)
        self.fig.canvas.get_tk_widget().grid(row = 1, column = 0, sticky = tk.NE + tk.SW)
        pinfo_label.grid(row = 2, column = 0, sticky = tk.E + tk.W)
        cframe.grid(row = 1, column = 1, rowspan = 2, sticky = tk.N + tk.S)
        #self.fig.canvas.get_tk_widget().pack(side = tk.LEFT, fill = tk.BOTH, expand = 1)
        #pinfo_label.pack(side = tk.BOTTOM, fill = tk.X)
        #cframe.pack(fill = tk.Y, expand = 1)

        print(uptime.uptime(), "Loading presets...")
        #FIXME: The first time preset is loaded, tmax, limits etc are not being
        # updated from some reason.
        self.controls['system']._load_preset('User defined')
        #self.controls['system']._load_preset('Lorentz attractor')

        self._init_keybindings()
        self.show_about()

    def _init_options(self, fname = None):
        opts = Options()
        opts.tmax           = tk.DoubleVar(self.root, 25)
        opts.dt             = tk.DoubleVar(self.root, 0.05)
        opts.reverse        = tk.BooleanVar(self.root, True)
        opts.quiver         = tk.BooleanVar(self.root, False)
        opts.nullclines     = tk.BooleanVar(self.root, False)
        opts.fixed_points   = tk.BooleanVar(self.root, False)
        opts.temporal       = tk.BooleanVar(self.root, False)
        opts.projection     = tk.StringVar(self.root, '2D')
        opts.limits         = Options()
        opts.limits.factor  = tk.DoubleVar(self.root, 2.5)
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

    def _set_limits(self):
        preset_name = self.controls['system'].preset.get()
        preset = presets.systems[preset_name]
        opts = self.opts
        x1, x2 = -5.0, 5.
        y1, y2 = -5.0, 5.
        z1, z2 = -5.0, 5.
        if preset is not None:
            if 'xlim' in preset:
                x1, x2 = preset['xlim']
            if 'ylim' in preset:
                y1, y2 = preset['ylim']
            if 'zlim' in preset:
                z1, z2 = preset['zlim']
        opts.limits.xmin.set(x1)
        opts.limits.xmax.set(x2)
        opts.limits.ymin.set(y1)
        opts.limits.ymax.set(y2)
        opts.limits.zmin.set(z1)
        opts.limits.zmax.set(z2)
        self._update_system_limits()

    def _handle_resize(self, *args):
        # This fixes a bug in 3d mode in which a projection of the trajectory
        # is drawn on the z=0 plane when window is resized.
        self.root.after_idle(self.update_fig)

    def _load_preset(self, name):
        opts = self.opts
        preset = presets.systems[name]
        defaults = dict(tmax = 25, dt = 0.05, projection = '2D', reverse = True)
        for opt in ('tmax', 'dt', 'projection', 'reverse'):
            if opt in preset:
                opts[opt].set(preset[opt])
            else:
                opts[opt].set(defaults[opt])
        self._set_limits()
        self._reset_fig()
        self._set_proj()
        if 'locations' in preset:
            for pos in preset['locations']:
                self.add_location(pos)
        #self.update_trajectories()

    def _update_system_limits(self, evt = None, prompt = False):
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
        #self.fig.draw(force = True)
        self.update_fig()
        for i in range(3):
            if periodic[i]:
                self.system.limits[i] = plot_limits[i]
            else:
                self.system.limits[i] = helpers.scale_domain(plot_limits[i], factor)
        self.system.periodic = periodic
        if PROJECTIONS[self.opts.projection.get()] == 'polar':
            # Ignore xlim (theta) in polar mode.
            self.system.limits[0] = None
            self.system.periodic[0] = True

    def show_poincare_dialog(self, *args):
        if len(self.trajectories) == 0:
            sys.stderr.write("Add a trajectory first.\n")
            return
        self.stop_traj_animation()
        l = self.opts.limits
        x1, x2 = l.xmin.get(), l.xmax.get()
        y1, y2 = l.ymin.get(), l.ymax.get()
        z1, z2 = l.zmin.get(), l.zmax.get()
        w = PWindow(self.root, self.trajectories[self.last_loc],
                ((x1, x2), (y1, y2), (z1, z2)),
                geometry = self.root.winfo_geometry())

    def update_trajectories(self, *args):
        picked = list(self.trajectories.keys())
        if len(picked) == 0:
            preset_name = self.controls['system'].preset.get()
            preset = presets.systems[preset_name]
            if 'locations' in preset:
                picked = preset['locations']
        self.trajectories.clear()
        self.fixed_points.clear()
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
            self.fig.draw_nullclines(field, linestyles = "dashed")
        if self.opts.quiver.get():
            self.fig.draw_quiver(field_quiv, width = 0.001, headwidth = 5, scale = 50)
        if self.opts.fixed_points.get():
            self.fig.draw_fp(*self.fixed_points)
        #NOTE: Add other bg elements like FP, LC, here before draw and save.
        self.fig.draw(force = True)
        self.fig.save()
        for t in self.trajectories.values():
            self.fig.add_trajectory(t)
            #self.fig.draw_trajectory(t)
        self.fig.draw()

    def _reset_fig(self, *args):
        self.stop_traj_animation(update = False)
        self.trajectories.clear()
        self.fixed_points.clear()
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
        proj = PROJECTIONS[self.opts.projection.get()]
        if proj.lower() == '3d':
            self.controls['nullclines'].configure(state = tk.DISABLED)
            self.controls['quiver'].configure(state = tk.DISABLED)
        else:
            self.controls['nullclines'].configure(state = tk.NORMAL)
            self.controls['quiver'].configure(state = tk.NORMAL)
        xlims, ylims, zlims = self.controls['limits']
        ylims[1].name = 'ymax:'
        for l in self.controls['limits']:
            l[0].enable()
            l[1].enable()
        if proj.lower() == 'polar':
            ylims[1].name = 'rmax:'
            for l in self.controls['limits']:
                l[0].disable()
                l[1].disable()
            self.controls['limits'][1][1].enable()
        self.fig.set_proj(proj)
        self._update_system_limits()
        self.update_fig()

    def anim_update(self):
        tstep = self.anim_tstep.get()
        anim_time = tstep * self.opts.dt.get()
        if anim_time > self.opts.tmax.get():
            anim_time = 0.0
            self.anim_tstep.set(0)
        else:
            self.anim_tstep.set(tstep+1)
        self.anim_info.set("t = %0.3f" % anim_time)
        self.pointer_info.set("t = %0.3f" % anim_time)
        self.fig.anim_update(anim_time, self.trajectories.values())

    def stop_traj_animation(self, update = True):
        self.anim_tstep.set(0)
        self.anim_timer.stop()
        self.anim_running.set(False)
        if not self.anim_running.get() and update:
            self.update_fig()
        self.controls['anim'].configure(text = 'Animate')
        self.anim_info.set("t = 0.000")

    def toggle_traj_animation(self):
        if len(self.trajectories) == 0:
            return
        if self.anim_running.get():
            print("Stopping animation...")
            self.anim_timer.stop()
            self.anim_running.set(False)
            self.controls['anim'].configure(text = 'Animate')
        else:
            print("Starting animation...")
            self.anim_timer.start()
            self.anim_running.set(True)
            self.controls['anim'].configure(text = 'Pause')

    def add_location(self, pos):
        styles = ["b", "g", "r", "c", "m", "y", "k"]
        threshold = 1e-4
        if len(pos) == 2:
            pos = pos[0], pos[1], 0.0
        pos = tuple(pos)
        #Search for a fixed point
        fp = self.system.find_fp(pos, threshold = 1e-4)
        if fp is not None:
            fp_clean = tuple(np.round(fp, 3))
            if not fp_clean in self.fixed_points:
                self.fixed_points.add(fp_clean)
                self.fig.draw_fp(fp_clean)
                self.update_fig()
                vel = self.system(fp_clean)
                print("Found a fixed point:\n", fp_clean, vel)
        if pos in self.trajectories:
            print("Trajectory already exists.")
            return
        try:
            t = Timer()
            t.start()
            traj = self.system.trajectory(pos, self.opts.tmax.get(), threshold = threshold,
                bidirectional = self.opts.reverse.get(), nsteps = 5 * (self.opts.tmax.get()/self.opts.dt.get()),
                max_step = self.opts.dt.get(), use_ode = True)
            t.stop()
            print("Computing trajectory (%d points) took %g seconds" % (len(traj.x), t.seconds()))
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
        traj.style = styles[style]
        self.fig.add_trajectory(traj)
        #self.fig.draw_trajectory(traj)
        self.fig.draw()
        self.last_loc = pos

    def button_press(self, evt):
        self.press_loc = evt.x, evt.y

    def mouse_move(self, evt):
        s = ""
        if evt.inaxes is None:
            s = ""
        elif evt.inaxes == self.fig.ax_polar:
            s = evt.inaxes.format_coord(evt.xdata, evt.ydata)
        elif evt.inaxes == self.fig.ax_rect:
            s = "x = %0.3f, y = %0.3f" % (evt.xdata, evt.ydata)
        elif evt.inaxes == self.fig.ax_3d:
            s = evt.inaxes.format_coord(evt.xdata, evt.ydata)
            if evt.button is None:
                x, y, z = [float(c.split('=')[-1]) for c in s.split(',')]
                s = "x = %0.3f, y = %0.3f, z = %0.3f" % (x, y, z)
        self.pointer_info.set(s)

    def button_release(self, evt):
        if not evt.inaxes is self.fig.ax_main:
            return
        if evt.button == 1:
            if self.fig.ax_main is self.fig.ax_3d:
                # Rotation or zoom
                self.update_fig()
        if self.anim_running.get():
            return
        release_loc = evt.x, evt.y
        d = helpers.distance(self.press_loc, release_loc)
        if d > 5:
            return
        x, y, z = evt.xdata, evt.ydata, 0.0
        if self.fig.ax_main is self.fig.ax_3d:
            s = evt.inaxes.format_coord(evt.xdata, evt.ydata)
            x, y, z = [float(c.split('=')[-1]) for c in s.split(',')]
        pos = x, y, z
        #print("Position:", pos)
        #print("Limits:", self.fig.get_limits())
        #if not helpers.is_inside(pos, self.fig.get_limits()):
        #    print("Position outside limits.")
        #    return
        self.add_location(pos)

    def _add_widgets(self, frame):
        controls = dict()
        f_system = DSFrame(frame, self.system, preset_cmd = self._load_preset,
                command = self.update_trajectories)
        f_system.grid(sticky = tk.W + tk.E)
        controls['system'] = f_system

        # Plot options frame
        f_controls = tk.LabelFrame(frame, text = "Options")
        f_controls.grid(sticky = tk.E + tk.W)
        row = 0
        f = tk.Frame(f_controls)
        f.grid(row=row, columnspan=2)
        tk.Label(f, text = "Projection:",
                anchor = tk.E).grid(row=0, column = 0, sticky = tk.E)
        optmenu = tk.OptionMenu(f, self.opts.projection,
                *PROJECTIONS.keys(), command = self._set_proj)
        optmenu.configure(width = 12)
        optmenu.grid(row = 0, column = 1, sticky = tk.W+tk.E)
        controls['proj'] = optmenu

        row += 1
        btn_nullc = tk.Checkbutton(f_controls, text = "Nullclines",
                variable = self.opts.nullclines, command = self.update_fig)
        btn_nullc.grid(row=row, column=0, sticky = tk.W)
        btn_quiver = tk.Checkbutton(f_controls, text = "Quiver",
                variable = self.opts.quiver, command = self.update_fig)
        btn_quiver.grid(row = row, column = 1, sticky = tk.W)
        controls['nullclines'] = btn_nullc
        controls['quiver'] = btn_quiver
        row += 1
        btn_temporal = tk.Checkbutton(f_controls, text = "Time Series",
                variable = self.opts.temporal, command = self._set_temporal)
        btn_temporal.grid(row=row, column=0, sticky = tk.W)
        btn_fp = tk.Checkbutton(f_controls, text = "Fixed Points",
                variable = self.opts.fixed_points, command = self.update_fig)
        btn_fp.grid(row = row, column = 1, sticky = tk.W)
        controls['temporal'] = btn_temporal
        controls['fp'] = btn_fp

        row += 1
        PEntry(f_controls, "Tstep:", self.opts.dt).grid(row = row, column = 0)
        PEntry(f_controls, 'Tmax:', self.opts.tmax).grid(row = row, column = 1)

        row += 1
        limits = self.opts.limits
        labels = [('xmin:', 'xmax:'), ('ymin:', 'ymax:'), ('zmin:', 'zmax:')]
        vars = [(limits.xmin, limits.xmax), (limits.ymin, limits.ymax), (limits.zmin, limits.zmax)]
        controls['limits'] = [[None, None],[None, None],[None, None]]
        for r in 0, 1, 2:
            for c in 0, 1:
                e = PEntry(f_controls, labels[r][c], vars[r][c], command = self._update_system_limits)
                e.grid(row = row+r, column = c)
                controls['limits'][r][c] = e
        row += 3
        tk.Button(f_controls, text = "Apply limits",
                command = self._update_system_limits).grid(row = row, column = 0)
        tk.Button(f_controls, text = 'Defaults',
                command = self._set_limits).grid(row = row, column = 1)

        f_anim = tk.LabelFrame(frame, text = 'Animation')
        f_anim.columnconfigure(0, weight=1)
        f_anim.grid(sticky = tk.E + tk.W)
        b_toggle = tk.Button(f_anim, text = "Animate", command = self.toggle_traj_animation,
                background = "#0000aa", activebackground = "#3333ff",
                foreground = "white", activeforeground = "white", font = "sans 16 bold",
                height = 1, width = 6)
        b_toggle.grid(row = 0, rowspan = 2, column = 0)
        b_stop = tk.Button(f_anim, text = "Stop", command = self.stop_traj_animation, font = "sans 10 bold",
                background = "#aa0000", activebackground = "#ff5555",
                foreground = "white", activeforeground = "white")
        b_stop.grid(row = 0, rowspan=2, column = 1)
        controls['anim'] = b_toggle

        f_update = tk.Frame(frame)
        f_update.columnconfigure(0, weight=1)
        f_update.grid(sticky = tk.E + tk.W + tk.S)
        # Allocate remaining vertical space to this frame.
        # FIXME: Can this be done better?
        frame.rowconfigure(f_update.grid_info()['row'], weight = 1)
        f_update.rowconfigure(0, weight=1)
        tk.Button(f_system, text = "Update", command = self.update_trajectories,
                background = "#0000aa", activebackground = "#3333ff",
                foreground = "white", activeforeground = "white", font = "sans 16 bold",
                height = 1, width = 6).grid(
                        columnspan=2, sticky = tk.S+tk.W+tk.E)
        tk.Button(f_update, text = "Clear", command = self._reset_fig, font = "sans 10 bold",
                background = "#aa0000", activebackground = "#ff5555",
                foreground = "white", activeforeground = "white").grid(
                        row = 0, columnspan = 2, sticky = tk.S)
        return controls

    def show_about(self):
        print("About this application")
        AboutDialog(self.root)

    def _add_menubar(self):
        menubar = tk.Menu(self.root)
        filemenu = tk.Menu(menubar, tearoff = False)
        menubar.add_cascade(label = 'File', menu=filemenu)
        filemenu.add_command(label = 'Print', command = self.save)
        filemenu.add_command(label = 'Quit', command = self.root.quit)
        viewmenu = tk.Menu(menubar, tearoff = False)
        menubar.add_cascade(label = 'View', menu=viewmenu)
        #viewmenu.add_checkbutton(label = 'Time series',
        #        variable = self.opts.temporal,
        #        command = self._set_temporal)
        #viewmenu.add_checkbutton(label = 'Nullclines',
        #        variable = self.opts.nullclines,
        #        command = self.update_fig)
        #viewmenu.add_checkbutton(label = 'Quiver',
        #        variable = self.opts.quiver,
        #        command = self.update_fig)
        #viewmenu.add_checkbutton(label = 'Fixed points',
        #        variable = self.opts.fixed_points,
        #        command = self.update_fig)
        #viewmenu.add_separator()
        #viewmenu.add_command(label = "Fullscreen", command = )
        #viewmenu.add_separator()
        viewmenu.add_command(label = "Poincare", command = self.show_poincare_dialog)
        helpmenu = tk.Menu(menubar, tearoff = False)
        menubar.add_cascade(label = 'Help', menu=helpmenu)
        helpmenu.add_command(label = 'About', command = self.show_about)
        return menubar

    def save(self):
        f = asksaveasfilename(defaultextension = ".pdf",
                parent = self.root, title = "Save as", initialfile = 'figure',
                filetypes = [("PDF files", "*.pdf")])
        f = str(f)
        if not f.endswith('.pdf'):
            return
        print("Saving to", f, type(f))
        self.fig.fig.savefig(f)

    def _init_keybindings(self):
        self.root.bind_all('<Control-KeyPress-p>', lambda *args: self.save())
        self.root.bind_all('<Control-KeyPress-s>', lambda *args: self.save())
        self.root.bind_all('<Control-KeyPress-q>', lambda *args: self.root.quit())


if __name__ == '__main__':
    np.seterr(invalid = 'print')
