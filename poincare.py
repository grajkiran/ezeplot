#!/usr/bin/python
try:
    import tkinter as tk
except:
    import Tkinter as tk
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from widgets import VEntry
import numpy as np

def sign(x):
    if x > 0: return 1
    if x < 0: return -1
    return 0

class Plane:
    """A plane in cartesian 3D space"""
    def __init__(self, a, b, c, d):
        self.a = float(a)
        self.b = float(b)
        self.c = float(c)
        self.d = float(d)

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
            #print(s_prev, s, p)
            if s_prev != s:
                #print("Intersection within:", p_prev, p)
                p_int = self.intersection(p_prev, p)
                if s_prev > 0:
                     top2bot.append(p_int)
                else:
                     bot2top.append(p_int)
            p_prev = p
            s_prev = s
        return top2bot, bot2top

class PWindow(tk.Toplevel):
    def __init__(self, master, trajectories):
        tk.Toplevel.__init__(self, master)
        self.trajectories = trajectories
        self.fig = matplotlib.figure.Figure()
        self.canvas = FigureCanvasTkAgg(self.fig, master = self)
        self.canvas.get_tk_widget().pack(side = tk.BOTTOM, fill = tk.BOTH, expand = 1)
        self.fig.subplots_adjust(left = 0.0, right = 1.0, bottom = 0.0, top = 1.0)
        self.ax = self.fig.add_subplot(111, projection = '3d')
        #self.ax.axis('off')
        self.a = tk.DoubleVar(self, 0.0)
        self.b = tk.DoubleVar(self, 1.0)
        self.c = tk.DoubleVar(self, 0.0)
        self.d = tk.DoubleVar(self, 0.0)
        cframe = tk.Frame(self)
        cframe.pack()
        self._add_widgets(cframe)

    def _add_widgets(self, frame):
        tk.Label(frame, text = "Equation of plane:").pack(side = tk.LEFT)
        VEntry(frame, textvariable = self.a, width = 5).pack(side = tk.LEFT)
        tk.Label(frame, text = "x + ").pack(side = tk.LEFT)
        VEntry(frame, textvariable = self.b, width = 5).pack(side = tk.LEFT)
        tk.Label(frame, text = "y + ").pack(side = tk.LEFT)
        VEntry(frame, textvariable = self.c, width = 5).pack(side = tk.LEFT)
        tk.Label(frame, text = "z + ").pack(side = tk.LEFT)
        VEntry(frame, textvariable = self.d, width = 5).pack(side = tk.LEFT)
        tk.Label(frame, text = " = 0").pack(side = tk.LEFT)

def toroidal_trajectory(R = 5.0, r = 3.0, w1 = 1.0, w2 = 10.0, tmax = 10.0, intervals = 1000):
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
    print "Top to bottom:"
    for p in top2bot:
        print p
    print "Bottom to top:"
    for p in bot2top:
        print p
    import matplotlib.pyplot as plt
    plt.plot(t.x, t.y, 'ro-')
    plt.show()

if __name__ == '__main__':
    #test_crossings()
    root = tk.Tk()
    w = PWindow(root, None)
    t = toroidal_trajectory()
    p = PSection(0.0, 2.0, 1.0, 0.0)
    down, up = p.compute_crossings(t)
    up = np.array(up)
    down = np.array(down)
    #w.ax.plot(t.x, t.y, t.z, 'r-')
    #w.ax.plot(up[:,0], up[:,1], up[:,2], 'b*')
    #w.ax.plot(down[:,0], down[:,1], down[:,2], 'go')
    th1, th2, x, y, z = torus()
    w.ax.plot_surface(x, y, z, rstride = 1, cstride = 1, shade = False)
    w.mainloop()
