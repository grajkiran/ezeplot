#!/usr/bin/python
try:
    import tkinter as tk
except:
    import Tkinter as tk
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
from widgets import VEntry
import numpy as np
from textwrap import wrap
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
        #print("Intersection:", dx, dy, dz)
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
                val += "%0.3g%s" % (abs(coeff), var)
        val += " = %0.3g" % -d
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
    def __init__(self, master, trajectory, limits, geometry = None):
        tk.Toplevel.__init__(self, master)
        if geometry is not None:
            self.geometry(geometry)
        self.transient(master)
        self.title("EzePlot - Poincare section")
        self.trajectory = trajectory
        self.fig = matplotlib.figure.Figure()
        self.canvas = FigureCanvasTkAgg(self.fig, master = self)
        #toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        #toolbar.update()
        self.canvas.get_tk_widget().pack(side = tk.LEFT, fill = tk.BOTH, expand = 1)
        self.fig.subplots_adjust(left = 0.0, right = 1.0, bottom = 0.0, top = 1.0)
        self.ax = Axes3D(self.fig)
        self.limits = limits
        self.p0_str = tk.StringVar(self, "0 0 0")
        self.p1_str = tk.StringVar(self, "0 0 1")
        self.factor = tk.DoubleVar(self, 0.5)
        self.draw_traj = tk.BooleanVar(self, True)
        self.auto_view = tk.BooleanVar(self, False)
        self.draw_plane = tk.BooleanVar(self, False)
        self.draw_endpoints = tk.BooleanVar(self, False)
        self.plane_eqn = tk.StringVar(self, "0.0 x + 0.0 y + z = 0")
        self.elev = tk.DoubleVar(self, 0.0)
        self.azim = tk.DoubleVar(self, 0.0)
        cframe = tk.Frame(self)
        cframe.pack()
        self._add_widgets(cframe)
        self.protocol('WM_DELETE_WINDOW', self.destroy)
        self.grab_set()
        self._update()
        self._plane_preset(2)
        self.wait_window(self)

    def _add_widgets(self, frame):
        f_plane = tk.LabelFrame(frame, text = "Plane")
        f_plane.grid(sticky = tk.E + tk.W)
        f_plane.columnconfigure(1, weight = 1)
        row = 0
        tk.Label(f_plane, text = 'Normal direction:').grid(row = row, columnspan = 2, sticky = tk.W)
        row += 1
        f_buttons = tk.Frame(f_plane)
        f_buttons.grid(row = row, columnspan = 2)#, sticky = tk.W+tk.E)
        tk.Button(f_buttons, text = "X", command = partial(self._plane_preset, 0)).pack(side = tk.LEFT)
        tk.Button(f_buttons, text = "Y", command = partial(self._plane_preset, 1)).pack(side = tk.LEFT)
        tk.Button(f_buttons, text = "Z", command = partial(self._plane_preset, 2)).pack(side = tk.LEFT)
        #tk.Button(f_buttons, text = "Screen", command = self._plane_preset).pack(side = tk.LEFT)
        #row += 1
        #tk.Label(f_plane, text = 'Two points:').grid(columnspan = 2, sticky = tk.W)
        row += 1
        tk.Label(f_plane, text = "From:").grid(row = row, column = 0,sticky = tk.E)
        VEntry(f_plane, self.p0_str, validator = parse_point,
                command = self._update).grid(row = row, column = 1, sticky = tk.E + tk.W)
        row += 1
        tk.Label(f_plane, text = "To:").grid(row = row, column = 0, sticky = tk.E)
        VEntry(f_plane, self.p1_str, validator = parse_point, width = 10,
                command = self._update).grid(row = row, column = 1, sticky = tk.E + tk.W)
        row += 1
        tk.Label(f_plane, text = "Ratio:").grid(row = row, column = 0, sticky = tk.S)
        tk.Scale(f_plane, orient = tk.HORIZONTAL, variable = self.factor,
                resolution = 0.005, from_ = 0.0, to = 1.0,
                command = self._update).grid(row = row, column= 1,
                        sticky = tk.E + tk.W)
        row += 1
        tk.Label(f_plane, text = "Equation of the plane:").grid(row = row, columnspan = 2, sticky = tk.W)
        row += 1
        tk.Label(f_plane, textvariable = self.plane_eqn, width = 30, anchor = tk.W,
                font = 'TkFixedFont').grid(row = row, columnspan = 2, sticky = tk.W)

        f_options = tk.LabelFrame(frame, text = "Options")
        f_options.grid(sticky = tk.E + tk.W)
        tk.Checkbutton(f_options, text = "Show trajectory", command = self._update,
                variable = self.draw_traj).grid(columnspan=2, sticky = tk.W)
        tk.Checkbutton(f_options, text = "Show end points", command = self._update,
                variable = self.draw_endpoints).grid(columnspan=2, sticky = tk.W)
        tk.Checkbutton(f_options, text = "Show plane", command = self._update,
                variable = self.draw_plane).grid(columnspan=2, sticky = tk.W)
        tk.Checkbutton(f_options, text = "Auto align", command = self._update,
                variable = self.auto_view).grid(columnspan=2, sticky = tk.W)

    def _plane_preset(self, direction = None):
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

    def _update(self, *args):
        t = self.trajectory
        p0 = parse_point(self.p0_str.get())
        p1 = parse_point(self.p1_str.get())
        plane = PSection.from_points(p0, p1, self.factor.get())
        self.plane_eqn.set(str(plane))
        if helpers.distance(p0, p1) == 0.0:
            return
        from_top, from_bot = plane.compute_crossings(t)
        xlim, ylim, zlim = self.limits
        self.ax.clear()
        #self.ax.set_axis_off()
        self.ax.set_frame_on(False)
        self.ax.grid(False)
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
            plane.plot(self.ax, self.limits, color = 'green',
                    shade = False, linewidth = 0.1, edgecolor = 'green')
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
    w = PWindow(root, t, ((-7,10), (-9,12), (-8,11)))
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
