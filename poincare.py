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
    from tkinter.filedialog import asksaveasfilename
except:
    import Tkinter as tk
    from tkFileDialog import asksaveasfilename
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
from widgets import VEntry
import numpy as np
from textwrap import wrap
import logging
from functools import partial
import helpers

def sign(x):
    if x > 0: return 1
    if x < 0: return -1
    return 0

def parse_point(s):
    xs, ys, zs = s.replace(',', ' ').split()
    return float(xs), float(ys), float(zs)
def point_to_str(p):
    p_strs = ["%0.3g" % x if np.abs(x) > .001 else "0" for x in p]
    return ", ".join(p_strs)

def elev_azim(p0, p1):
    dx = p1[0] - p0[0]
    dy = p1[1] - p0[1]
    dz = p1[2] - p0[2]
    r = np.sqrt(dx*dx + dy*dy + dz*dz)
    elev_r = np.arcsin(dz/r)
    azim_r = np.arctan2(dy, dx)
    elev_d = np.degrees(elev_r)
    azim_d = np.degrees(azim_r)
    return elev_d, azim_d

class Plane:
    """A plane in cartesian 3D space"""
    def __init__(self, a, b, c, d):
        self.a = float(a)
        self.b = float(b)
        self.c = float(c)
        self.d = float(d)
        self.coeffs = a, b, c, d

    @classmethod
    def from_points(cls, p0, p1, factor = 0.0):
        x0, y0, z0 = p0
        x1, y1, z1 = p1
        a = x1 - x0
        b = y1 - y0
        c = z1 - z0
        f = -factor
        f_1 = factor - 1
        d  = a*(f_1 * x0 + f * x1)
        d += b*(f_1 * y0 + f * y1)
        d += c*(f_1 * z0 + f * z1)
        return cls(a, b, c, d)

    def dot(self, p):
        x, y, z = map(float, p)
        return self.a*x + self.b*y + self.c*z + self.d

    def intersection(self, p1, p2):
        """Returns the coordinates of intersection of the line formed by points
        p1 and p2 with the plane."""
        x1, y1, z1 = map(float, p1)
        x2, y2, z2 = map(float, p2)
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1
        a, b, c, d = self.a, self.b, self.c, self.d
        den = (a*dx + b*dy + c*dz)
        if den == 0:
            raise RuntimeError("Could not compute intersection. Line appears to be parallel to the plane.")
        t = - (a*x1 + b*y1 + c*z1 + d)/den
        if t < 0 or t > 1:
            raise RuntimeError("Intersection does not fall within.")
        x = x1 + t * dx
        y = y1 + t * dy
        z = z1 + t * dz
        return (x, y, z)
    def plot(self, ax, limits, **kwargs):
        xlim, ylim, zlim = limits
        p = max(self.a, self.b, self.c)
        if p == 0.0:
            return
        x = np.linspace(xlim[0], xlim[1], 12)
        y = np.linspace(ylim[0], ylim[1], 12)
        z = np.linspace(zlim[0], zlim[1], 12)
        if p == self.a: # plot x vs y, z
            y, z = np.meshgrid(y, z)
            x = (-self.d - self.b*y - self.c*z)/self.a
        elif p == self.b: # plot y vs x, z
            x, z = np.meshgrid(x, z)
            y = (-self.d - self.a*x - self.c*z)/self.b
        elif p == self.c:
            x, y = np.meshgrid(x, y)
            z = (-self.d - self.a*x - self.b*y)/self.c
        else: # Should not come here.
            return
        return ax.plot_surface(x, y, z, **kwargs)

    def __str__(self):
        val = ""
        a, b, c, d = self.coeffs
        for coeff, var in (a, 'X'), (b, 'Y'), (c, 'Z'):
            if coeff == 0:
                continue
            if coeff < 0:
                val += " - "
            else:
                val += " + "
            if abs(coeff) == 1.0:
                val += "%s" % var
            else:
                val += "%r%s" % (abs(coeff), var)
        val += " = %r" % -d
        return "\n\t".join(wrap(val.strip().lstrip('+'), 36))

class PSection(Plane):
    def compute_crossings(self, trajectory):
        top2bot = []
        bot2top = []
        x, y, z = trajectory.x, trajectory.y, trajectory.z
        N = len(x)
        p_prev = (x[0], y[0], z[0])
        s_prev = sign(self.dot(p_prev))
        for i in range(1, N):
            p = x[i], y[i], z[i]
            s = sign(self.dot(p))
            if s_prev != s:
                p_int = self.intersection(p_prev, p)
                if s_prev > 0:
                     top2bot.append(p_int)
                else:
                     bot2top.append(p_int)
            p_prev = p
            s_prev = s
        return np.array(top2bot), np.array(bot2top)

class PWindow(tk.Toplevel):
    def __init__(self, app):
        tk.Toplevel.__init__(self, app.root)
        if app.icon is not None:
            self.tk.call("wm", "iconphoto", self._w, app.icon)
        self.app = app
        self.geometry(app.root.geometry())
        self.transient(app.root)
        self.title("EzePlot - Poincare section")
        self.trajectory = app.trajectories[app.last_loc]
        self.fig = matplotlib.figure.Figure()
        self.canvas = FigureCanvasTkAgg(self.fig, master = self)
        #toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        #toolbar.update()
        self.canvas.get_tk_widget().pack(side = tk.LEFT, fill = tk.BOTH, expand = 1)
        #self.fig.subplots_adjust(left = 0.0, right = 1.0, bottom = 0.0, top = 1.0)
        self.ax = self.fig.add_subplot(111, projection = '3d')
        self.ax_x = self.fig.add_subplot(222, label = 'x', visible = False)
        self.ax_y = self.fig.add_subplot(223, label = 'y', visible = False)
        self.ax_z = self.fig.add_subplot(224, label = 'z', visible = False)
        self.mode = 1
        self.set_mode(1)
        self.limits = app._get_limits()
        self.p0_str = tk.StringVar(self, "0 0 0")
        self.p1_str = tk.StringVar(self, "0 0 1")
        self.factor = tk.DoubleVar(self, 0.5)
        self.plane_direction = tk.IntVar(self, 0)
        self.draw_traj = tk.BooleanVar(self, True)
        self.auto_view = tk.BooleanVar(self, False)
        self.draw_plane = tk.BooleanVar(self, True)
        self.draw_endpoints = tk.BooleanVar(self, False)
        self.first_returns = tk.BooleanVar(self, False)
        self.plane_eqn = tk.StringVar(self, "0.0 x + 0.0 y + z = 0")
        self.elev = tk.DoubleVar(self, 0.0)
        self.azim = tk.DoubleVar(self, 0.0)
        cframe = tk.Frame(self)
        cframe.pack(expand = True, fill = tk.Y)
        self.controls = self._add_widgets(cframe)
        self.menu = self._add_menubar()
        self.config(menu = self.menu)
        self.protocol('WM_DELETE_WINDOW', self.destroy)
        self.grab_set()
        self._update()
        self._plane_preset()
        self.save_section('poincare.txt')
        self.wait_window(self)

    def _add_widgets(self, frame):
        controls = dict()
        f_plane = tk.Frame(frame)#, text = "Plane")
        f_plane.grid(sticky = tk.E + tk.W)
        f_plane.columnconfigure(1, weight = 1)
        row = 0
        tk.Label(f_plane, text = 'Normal to plane',
                justify = tk.LEFT).grid(row = row, columnspan = 2, sticky = tk.W)
        row += 1
        f_buttons = tk.Frame(f_plane)
        f_buttons.grid(row = row, columnspan = 2)#, sticky = tk.W+tk.E)
        tk.Radiobutton(f_buttons, text = "X", variable = self.plane_direction,
                value = 0, command = self._plane_preset).pack(side = tk.LEFT)
        tk.Radiobutton(f_buttons, text = "Y", variable = self.plane_direction,
                value = 1, command = self._plane_preset).pack(side = tk.LEFT)
        tk.Radiobutton(f_buttons, text = "Z", variable = self.plane_direction,
                value = 2, command = self._plane_preset).pack(side = tk.LEFT)
        #tk.Button(f_buttons, text = "Screen", command = self._plane_preset).pack(side = tk.LEFT)
        #row += 1
        #tk.Label(f_plane, text = 'Two points:').grid(columnspan = 2, sticky = tk.W)
        #row += 1
        #tk.Label(f_plane, text = "From:").grid(row = row, column = 0,sticky = tk.E)
        #VEntry(f_plane, self.p0_str, validator = parse_point,
        #        command = self._update).grid(row = row, column = 1, sticky = tk.E + tk.W)
        #row += 1
        #tk.Label(f_plane, text = "To:").grid(row = row, column = 0, sticky = tk.E)
        #VEntry(f_plane, self.p1_str, validator = parse_point, width = 10,
        #        command = self._update).grid(row = row, column = 1, sticky = tk.E + tk.W)
        row += 1
        tk.Label(f_plane, text = "Location of plane").grid(row = row, column = 0, sticky = tk.S)
        scale = tk.Scale(f_plane, orient = tk.HORIZONTAL, variable = self.factor,
                resolution = 0.005, from_ = 0.0, to = 1.0,
                command = self._update)
        scale.grid(row = row+1, columnspan=2,
                        sticky = tk.E + tk.W)
        controls['location'] = scale
        #row += 1
        #tk.Label(f_plane, text = "Equation of the plane:").grid(row = row, columnspan = 2, sticky = tk.W)
        #row += 1
        #tk.Label(f_plane, textvariable = self.plane_eqn, width = 30, anchor = tk.W,
        #        font = 'TkFixedFont').grid(row = row, columnspan = 2, sticky = tk.W)

        f_options = tk.LabelFrame(frame, text = "Options")
        f_options.grid(sticky = tk.E + tk.W)
        tk.Checkbutton(f_options, text = "Show trajectory", command = self._update,
                variable = self.draw_traj).grid(columnspan=2, sticky = tk.W)
        #tk.Checkbutton(f_options, text = "Show end points", command = self._update,
        #        variable = self.draw_endpoints).grid(columnspan=2, sticky = tk.W)
        tk.Checkbutton(f_options, text = "Show plane", command = self._update,
                variable = self.draw_plane).grid(columnspan=2, sticky = tk.W)
        #tk.Checkbutton(f_options, text = "Auto align", command = self._update,
        #        variable = self.auto_view).grid(columnspan=2, sticky = tk.W)
        #tk.Checkbutton(f_options, text = "First return maps", command = self._update,
        #        variable = self.first_returns).grid(columnspan = 2, sticky = tk.W)
        tk.Label(frame, text = "Blue: Direction of\n         plane normal", fg = "blue",
                justify = tk.LEFT).grid()
        tk.Label(frame, text = "Red: Opposite to\n        plane normal", fg = "red",
                justify = tk.LEFT).grid()

        close_btn = tk.Button(frame, text = "Close", command = self.destroy, font = "sans 10 bold",
                background = "#aa0000", activebackground = "#ff5555",
                foreground = "white", activeforeground = "white")
        close_btn.grid(columnspan = 2, sticky = tk.S)
        frame.rowconfigure(close_btn.grid_info()['row'], weight = 1)

        return controls

    def _add_menubar(self):
        menubar = tk.Menu(self)
        filemenu = tk.Menu(menubar, tearoff = False)
        menubar.add_cascade(label = 'File', menu=filemenu)
        filemenu.add_command(label = 'Save Poincare section data', command = self.save_section)
        filemenu.add_command(label = 'Print', command = self.save)
        filemenu.add_separator()
        filemenu.add_command(label = 'Close', command = self.destroy)
        return menubar

    def save(self):
        f = asksaveasfilename(defaultextension = ".pdf",
                parent = self, title = "Save as", initialfile = 'figure',
                filetypes = [("PDF files", "*.pdf")])
        f = str(f)
        if not f.endswith('.pdf'):
            return
        logging.info("Saving to %s" % f)
        self.fig.savefig(f)

    def save_section(self, fname = None):
        if fname is None:
            f = asksaveasfilename(defaultextension = ".txt",
                    parent = self, title = "Save as", initialfile = 'poincare',
                    filetypes = [("Text files", "*.txt")])
            fname = str(f)
            if not fname.endswith('.txt'):
                return
        logging.info("Saving Poincare secion to %s" % fname)
        with open(fname, "w") as out:
            out.write("% Poincare section data\n")
            self.app.print_info(out)
            out.write("% Equation of the intersecting plane:\n")
            out.write("%%%30s\n" % (self.plane_eqn.get()))
            top, bottom = self.crossings
            out.write("% Points in the direction of plane normal\n")
            out.write("%% %-21s\t%-23s\t%-23s\n%%\n" % ("X", "Y", "Z"))
            for x, y, z in bottom:
                out.write("%-23r\t%-23r\t%-23r\n" % (x, y, z))
            out.write("%\n%\n% Points in the direction opposite to plane normal\n")
            out.write("%% %-21s\t%-23s\t%-23s\n%%\n" % ("X", "Y", "Z"))
            for x, y, z in top:
                out.write("%-23r\t%-23r\t%-23r\n" % (x, y, z))

    def _plane_preset(self):
        direction = self.plane_direction.get()
        limits = self.ax.get_xlim(), self.ax.get_ylim(), self.ax.get_zlim()
        center = [sum(limits[i])/2.0 for i in range(3)]
        scale = self.controls['location']
        scale.configure(from_ = limits[direction][0], to = limits[direction][1])
        self.factor.set(center[direction])
        self._update()

    def _plane_preset_old(self, direction = None):
        elevs = [0.0, 0.0, 90.0]
        azims = [0.0, 90.0, 0.0]
        limits = self.ax.get_xlim(), self.ax.get_ylim(), self.ax.get_zlim()
        if direction is None:
            #From screen
            elev = np.radians(self.ax.elev)
            azim = np.radians(self.ax.azim)
        else:
            elev = np.radians(elevs[direction])
            azim = np.radians(azims[direction])
        center = [sum(limits[i])/2.0 for i in range(3)]
        c1 = [l[0] for l in limits]
        c2 = [l[1] for l in limits]
        r = 0.5 * helpers.distance(c1, c2)
        dz = r * np.sin(elev)
        dy = r * np.sin(azim)*np.cos(elev)
        dx = r * np.cos(azim)*np.cos(elev)
        p0 = [center[0]-dx, center[1]-dy, center[2]-dz]
        p1 = [center[0]+dx, center[1]+dy, center[2]+dz]
        self.p0_str.set(point_to_str(p0))
        self.p1_str.set(point_to_str(p1))
        self.factor.set(0.5)
        self._update()

    def set_mode(self, mode = 1):
        if mode == 1:
            for ax in self.ax_x, self.ax_y, self.ax_z:
                ax.set_visible(False)
                self.fig.delaxes(ax)
            self.ax.change_geometry(1, 1, 1)
        elif mode == 4:
            for ax in self.ax_x, self.ax_y, self.ax_z:
                self.fig.add_axes(ax)
                ax.set_visible(True)
            self.ax.change_geometry(2, 2, 1)
        else:
            raise NotImplementedError("Unsupported mode: " + str(mode))
        self.mode = mode

    def _update_first_returns(self, from_top, from_bot):
        axes = (self.ax_x, self.ax_y, self.ax_z)
        labels = ('X', 'Y', 'Z')
        for i in range(3):
            axes[i].clear()
            axes[i].set_xlim(self.limits[i])
            axes[i].set_ylim(self.limits[i])
            axes[i].grid()
            axes[i].set_xlabel(r'$%s_i$' % labels[i])
            axes[i].set_ylabel(r'$%s_{i+1}$' % labels[i])
            if len(from_top) > 2:
                data_top = from_top[:,i]
                axes[i].plot(data_top[:-1], data_top[1:], 'r.')
            if len(from_bot) > 2:
                data_bot = from_bot[:,i]
                axes[i].plot(data_bot[:-1], data_bot[1:], 'b.')

    def _update(self, *args):
        if self.first_returns.get():
            if self.mode != 4:
                self.set_mode(4)
        elif self.mode != 1:
            self.set_mode(1)
        t = self.trajectory
        coeffs = [0.0, 0.0, 0.0, 0.0]
        coeffs[3] = -self.factor.get()
        coeffs[self.plane_direction.get()] = 1.0
        plane = PSection(*coeffs)
        self.plane_eqn.set(str(plane))
        self.crossings = from_top, from_bot = plane.compute_crossings(t)
        if self.first_returns.get():
            self._update_first_returns(from_top, from_bot)
        xlim, ylim, zlim = self.limits
        self.ax.clear()
        self.ax.text2D(0.5, 1.05, "Poincare section", size = 16, weight = 'bold',
                transform = self.ax.transAxes, ha = 'center')
        self.ax.set_frame_on(False)
        self.ax.grid(False)
        self.ax.mouse_init(zoom_btn = [], rotate_btn = [1])
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        self.ax.set_zlim(zlim)
        if self.auto_view.get():
            elev, azim = elev_azim(p0, p1)
            self.ax.elev = elev
            self.ax.azim = azim
        if self.draw_traj.get():
            self.ax.plot(t.x, t.y, t.z, 'k--')
        if len(from_top) > 0:
            self.ax.plot(from_top[:,0], from_top[:,1], from_top[:,2], 'r.')
        if len(from_bot) > 0:
            self.ax.plot(from_bot[:,0], from_bot[:,1], from_bot[:,2], 'b.')
        if self.draw_endpoints.get():
            self.ax.plot([p0[0]], [p0[1]], [p0[2]], 'bo')
            self.ax.plot([p1[0]], [p1[1]], [p1[2]], 'ro')
            self.ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]], 'k-')
        if self.draw_plane.get():
            plane.plot(self.ax, self.limits, color = 'yellow',
                    shade = False, linewidth = 0.1, edgecolor = 'yellow')
        self.canvas.draw()

def toroidal_trajectory(R = 5.0, r = 3.0, w1 = np.pi/3, w2 = 1.0,
        tmax = 1000.0, intervals = 5000):
    class traj: pass
    times = np.linspace(0, tmax, intervals)
    th1 = w1 * times
    th2 = w2 * times
    x = (R + r*np.cos(th2))*np.cos(th1)
    y = (R + r*np.cos(th2))*np.sin(th1)
    z = r * np.sin(th2)
    t = traj()
    traj.t = times
    traj.x = x
    traj.y = y
    traj.z = z
    return traj

def torus(R = 2, r = 1, intervals = 50):
    t1 = np.linspace(0, 2*np.pi, intervals)
    t2 = np.linspace(0, 2*np.pi, intervals)
    th1, th2 = np.meshgrid(t1, t2)
    x = (R + r*np.cos(th2))*np.cos(th1)
    y = (R + r*np.cos(th2))*np.sin(th1)
    z = r * np.sin(th2)
    print(len(th1))
    return th1, th2, x, y, z

def test_crossings():
    import math
    def traj_circle(radius, th_max = 10*math.pi, intervals = 100):
        import numpy as np
        import dynsystem
        th = np.linspace(0, th_max, intervals)
        s = radius * th
        x = radius * np.cos(th)
        y = radius * np.sin(th)
        z = np.zeros_like(th)
        return dynsystem.Trajectory(np.array([th, s, x, y, z]).transpose())
    p = PSection(2, 1, 4, 70)
    t = traj_circle(5.0, intervals = 5000)
    top2bot, bot2top = p.compute_crossings(t)
    print("Top to bottom:")
    for p in top2bot:
        print(p)
    print("Bottom to top:")
    for p in bot2top:
        print(p)
    import matplotlib.pyplot as plt
    plt.plot(t.x, t.y, 'ro-')
    plt.show()

if __name__ == '__main__':
    #test_crossings()
    root = tk.Tk()
    t = toroidal_trajectory()
    w = PWindow(root, t, ((-7,10), (-9,12), (-8,11)), geometry = "800x600")
    #p = PSection(0.0, 0.0, 1.0, 0.0)
    #down, up = p.compute_crossings(t)
    #up = np.array(up)
    #down = np.array(down)
    #w.ax.plot(t.x, t.y, t.z, 'r-')
    #w.ax.plot(up[:,0], up[:,1], up[:,2], 'b*')
    #w.ax.plot(down[:,0], down[:,1], down[:,2], 'go')
    #th1, th2, x, y, z = torus()
    #w.ax.plot_surface(x, y, z, rstride = 1, cstride = 1, shade = False)
    #w.mainloop()
