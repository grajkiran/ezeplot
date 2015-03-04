#!/usr/bin/python

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

POLAR2D = 0
POLAR3D = 1
CART3D = 2
CART2D = 3

class FigureWindow:
    def __init__(self, mode = 1):
        self.fig = plt.figure()
        self._blit = self.fig.canvas.supports_blit
        self.ax_main = self.fig.add_subplot(111, label = 'main', xlim = (-5, 5), ylim=(-5, 5))
        self.ax_x = self.fig.add_subplot(222, label = 'x', visible = False)
        self.ax_y = self.fig.add_subplot(223, label = 'y', visible = False, sharex = self.ax_x)
        self.ax_z = self.fig.add_subplot(224, label = 'z', visible = False, sharex = self.ax_x)
        self.bg = dict()
        # Dict of artists
        self.artists = dict(nullx = None, nully = None, quiver = None,
                picks = None, pquiver = None)
        self.trajectories = list()
        self.fig.show()
        self.set_mode(mode)

    def bind(self, event, handler):
        self.fig.canvas.mpl_connect(event, handler)

    def reset(self):
        self._discard_artists()
        # TODO: Empty the trajectories dict as well.
        self.update_bg()

    def save(self, axes):
        if self._blit:
            for ax in axes:
                self.bg[ax] = self.fig.canvas.copy_from_bbox(ax.bbox)

    def restore(self, axes):
        if self._blit:
            for ax in axes:
                if ax in self.bg:
                    self.fig.canvas.restore_region(self.bg[ax])

    def _discard_artists(self, *names):
        """Remove the specified artists from the list to be updated on each draw"""
        if len(names) == 0:
            names = self.artists.keys()
        for name in names:
            if name in self.artists and self.artists[name] is not None:
                self.artists[name].set_visible(False)
                self.artists[name].remove()
                self.artists[name] = None

    def _draw_artist(self, artist, ax = None):
        """Draws the given artists"""
        if artist is None:
            return
        if ax is None:
            ax = self.ax_main
        if not artist in ax.get_children():
            ax.add_artist(artist)
        if self._blit:
            ax.draw_artist(artist)

    def update(self):
        if self._blit:
            self.fig.canvas.blit()
        else:
            self.fig.canvas.draw_idle()

    def set_proj(self, projection = 'rectilinear'):
        geom = self.ax_main.get_geometry()
        self.fig.delaxes(self.ax_main)
        self.ax_main = self.fig.add_subplot(*geom, projection = projection)
        self.reset()

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
        self.update_bg()

    def update_bg(self, *axes):
        if len(axes) == 0:
            axes = self.fig.get_axes()
        for a in axes:
            a.clear()
            a.set_xlim(a.get_xlim())
            a.set_ylim(a.get_ylim())
            a.grid(True)
        for artist in self.artists.values():
            if artist is not None:
                self._draw_artist(artist)
        self.fig.canvas.draw()
        self.save(axes)

    def nullclines(self, field = None):
        """Plots the nullclines of the given field and stores the results.
        If field is None, then existing nullclines are removed. Otherwise,
        existing nullclines are updated.
        """
        self._discard_artists('nullx', 'nully')
        # Create nullclines and set self.artists['nullclines']
        if field is not None:
            x, y, u, v = field
            self.artists['nullx'] = self.ax_main.contour(x, y, u, [0.0]).collections[0]
            self.artists['nully'] = self.ax_main.contour(x, y, v, [0.0]).collections[0]
        self.update_bg()

    def quiver(self, field = None, unscaled = False, **opts):
        self._discard_artists('quiver')
        if field is not None:
            x, y, u, v = field
            if unscaled:
                mag = np.ma.masked_less(np.sqrt(u*u+v*v), 1.0e-4)
                u = u/mag
                v = v/mag
            self.artists['quiver'] = self.ax_main.quiver(x, y, u, v, **opts)
        self.update_bg()

    def add_trajectory(self, traj, style = 'r-', mfc = 'none', marker = 'o'):
        self.trajectories.append(traj)
        traj.line, = self.ax_main.plot(traj.x, traj.y, style)
        traj.marker_start, = self.ax_main.plot([traj.start[0]], [traj.start[1]], marker, mfc = mfc)
        self._draw_artist(traj.line)
        self._draw_artist(traj.marker_start)
        self.update()

def test_fw(root):
    import functools, itertools
    import dynsystem
    import numpy as np
    fw = FigureWindow()
    system = dynsystem.DynamicSystem('y', '-x + mu * (1 - x*x)*y', dict(mu = 1.0))
    x, y = np.meshgrid(np.linspace(-5, 5, 25), np.linspace(-5, 5, 25))
    u, v = system((x,y))
    field = x, y,  u, v
    mode = itertools.cycle([1, 4])
    proj = itertools.cycle(['rectilinear', 'polar', '3d'])
    v_nullc = tk.IntVar(root, 0)
    v_quiver = tk.IntVar(root, 0)
    v_quiv_unscaled = tk.IntVar(root, 1)
    next(mode), next(proj)
    def handle_click(evt):
        pos = evt.xdata, evt.ydata
        traj = system.trajectory(pos, dt = 0.05, max_steps = 500, direction = 0)
        fw.add_trajectory(traj)
    def toggle_mode():
        fw.set_mode(next(mode))
    def toggle_proj():
        fw.set_proj(next(proj))
    def h_nullc():
        if v_nullc.get():
            fw.nullclines(field)
        else:
            fw.nullclines()
    def h_quiver():
        if v_quiver.get():
            fw.quiver(field, unscaled = v_quiv_unscaled.get())
        else:
            fw.quiver()
    tk.Button(root, text = 'mode', command = toggle_mode).pack()
    tk.Button(root, text = 'projection', command = toggle_proj).pack()
    tk.Checkbutton(root, text = "Nullclines", variable = v_nullc, command = h_nullc).pack()
    tk.Checkbutton(root, text = "Quiver", variable = v_quiver, command = h_quiver).pack()
    tk.Checkbutton(root, text = "Unscaled", variable = v_quiv_unscaled).pack()
    fw.bind('button_release_event', handle_click)

if __name__ == '__main__':
    try:
        import tkinter as tk
    except:
        import Tkinter as tk
    root = tk.Tk()
    root.protocol('WM_DELETE_WINDOW', root.quit)
    test_fw(root)
    root.mainloop()
