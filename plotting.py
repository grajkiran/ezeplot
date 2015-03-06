#!/usr/bin/python

try:
    import tkinter as tk
except:
    import Tkinter as tk
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import helpers

matplotlib.rcParams['toolbar'] = 'None'

class Figure:
    def __init__(self, mode = 1, blit = True, master = None):
        if master is not None:
            self.fig = matplotlib.figure.Figure()
            self.canvas = FigureCanvasTkAgg(self.fig, master = master)
            self.canvas.get_tk_widget().pack(side = tk.LEFT, fill = tk.BOTH, expand = 1)
        else:
            self.fig = plt.figure()
            self.canvas = self.fig.canvas
        self._blit = blit and self.fig.canvas.supports_blit
        self.ax_main = self.fig.add_subplot(111, label = 'main', xlim = (-5, 5), ylim=(-5, 5))
        self.ax_x = self.fig.add_subplot(222, label = 'x', visible = False)
        self.ax_y = self.fig.add_subplot(223, label = 'y', visible = False, sharex = self.ax_x)
        self.ax_z = self.fig.add_subplot(224, label = 'z', visible = False, sharex = self.ax_x)
        self.bg = dict()
        self.trajectories = []
        if master is None:
            self.fig.show()
        self.set_mode(mode)

    def get_diagonal(self):
        lim = self.get_limits()
        p1, p2 = np.array(lim).transpose()
        return helpers.distance(p1, p2)

    def get_limits(self):
        xlim = self.ax_main.get_xlim()
        ylim = self.ax_main.get_ylim()
        if not hasattr(self.ax_main, 'get_zlim'):
            return xlim, ylim
        zlim = self.ax_main.get_zlim()
        return xlim, ylim, zlim

    def bind(self, event, handler):
        self.fig.canvas.mpl_connect(event, handler)

    def save(self, *axes):
        if not self._blit: return
        if len(axes) == 0:
            axes = self.fig.get_axes()
        for ax in axes:
            self.bg[ax] = self.fig.canvas.copy_from_bbox(ax.bbox)

    def restore(self, *axes):
        if not self._blit: return
        if len(axes) == 0:
            axes = self.fig.get_axes()
        for ax in axes:
            if ax in self.bg:
                self.fig.canvas.restore_region(self.bg[ax])

    def _draw_artist(self, artist, ax = None):
        """Draws the given artists. First adds it to the list of artists in the axes"""
        if artist is None:
            return
        if ax is None:
            ax = self.ax_main
        if not artist in ax.get_children():
            ax.add_artist(artist)
        if self._blit:
            ax.draw_artist(artist)

    def set_proj(self, projection = 'rectilinear'):
        geom = self.ax_main.get_geometry()
        self.fig.delaxes(self.ax_main)
        self.ax_main = self.fig.add_subplot(*geom, projection = projection)
        self.draw(force = True)

    def set_mode(self, mode = 1):
        if mode == 1:
            for ax in self.ax_x, self.ax_y, self.ax_z:
                ax.set_visible(False)
                self.fig.delaxes(ax)
            self.ax_main.change_geometry(1, 1, 1)
        elif mode == 4:
            for ax in self.ax_x, self.ax_y, self.ax_z:
                self.fig.add_axes(ax)
                ax.set_visible(True)
            self.ax_main.change_geometry(2, 2, 1)
        else:
            raise NotImplementedError("Unsupported mode: " + str(mode))

    def clear(self):
        for a in self.fig.get_axes():
            a.clear()
            a.set_xlim(a.get_xlim())
            a.set_ylim(a.get_ylim())
            a.grid(True)

    def draw(self, force = False):
        if force:
            self.fig.canvas.draw()
        if self._blit:
            self.fig.canvas.blit()
        else:
            self.fig.canvas.draw_idle()

    # NOTE: the draw_* functions must not call draw or blit. 
    def draw_nullclines(self, field, **kwargs):
        x, y, u, v = field
        self.ax_main.contour(x, y, u, [0.0], **kwargs)
        self.ax_main.contour(x, y, v, [0.0], **kwargs)

    def draw_quiver(self, field, scaled = True, **kwargs):
        x, y, u, v = field
        if scaled:
            mag = np.ma.masked_less(np.sqrt(u*u+v*v), 1.0e-4)
            u = u/mag
            v = v/mag
        self.ax_main.quiver(x, y, u, v, **kwargs)

    def draw_trajectory(self, traj, t_anim = None):
        if t_anim is None:
            traj.marker_start.set_visible(True)
            traj.marker_anim.set_visible(False)
            traj.arrow.set_visible(True)
            traj.line.set_xdata(traj.x)
            traj.line.set_ydata(traj.y)
        else: # During animation
            if t_anim > traj.t[-1]:
                t_anim = traj.t[-1]
            last = traj.prev_index(t_anim)+1
            traj.line.set_xdata(traj.x[:last])
            traj.line.set_ydata(traj.y[:last])
            x = traj.x_ip(t_anim)
            y = traj.y_ip(t_anim)
            traj.marker_anim.set_xdata([x])
            traj.marker_anim.set_ydata([y])
            traj.marker_anim.set_visible(True)
        self._draw_artist(traj.line)
        self._draw_artist(traj.arrow)
        self._draw_artist(traj.marker_start)
        self._draw_artist(traj.marker_anim)

    def add_trajectory(self, traj, *args, **kwargs):
        self.trajectories.append(traj)
        self.init_trajectory(traj, *args, **kwargs)

    def init_trajectory(self, traj, style = 'r-', mfc = 'none', marker = 'o'):
        """Creates the line, start marker, arrow and animation marker"""
        traj.line, = self.ax_main.plot(traj.x, traj.y, style)
        traj.marker_start, = self.ax_main.plot([traj.start[0]], [traj.start[1]], marker, mfc = mfc)
        traj.marker_anim, = self.ax_main.plot([traj.start[0]], [traj.start[1]], marker)
        traj.marker_anim.set_visible(False)
        # Add arrow at about 10% diagonal length from the start location
        diag = self.get_diagonal()
        arrow_len = diag / 1000.0
        arrow_start = min(self.get_diagonal()/10.0, traj.dist[-1]/4.0)
        arrow_end = arrow_start + arrow_len
        t_start, t_end = traj.t_ip([arrow_start, arrow_end])
        x1, x2 = traj.x_ip([t_start, t_end])
        y1, y2 = traj.y_ip([t_start, t_end])
        traj.arrow = self.ax_main.annotate("", xy = (x2, y2), xytext = (x1, y1),
                arrowprops = dict(arrowstyle = "->"))

    def anim_update(self, time):
        self.restore(self.ax_main)
        for traj in self.trajectories:
            self.draw_trajectory(traj, time)
        self.draw()

if __name__ == '__main__':
    try:
        import tkinter as tk
    except:
        import Tkinter as tk
    root = tk.Tk()
    root.protocol('WM_DELETE_WINDOW', root.quit)
    test_fw(root)
    root.mainloop()
