#!/usr/bin/python

try:
    import tkinter as tk
except:
    import Tkinter as tk
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import numpy as np
import helpers

#matplotlib.rcParams['toolbar'] = 'None'

class Figure:
    def __init__(self, mode = 1, blit = True, master = None):
        if master is not None:
            self.fig = matplotlib.figure.Figure()
            self.canvas = FigureCanvasTkAgg(self.fig, master = master)
            self.canvas.get_tk_widget().pack(side = tk.LEFT, fill = tk.BOTH, expand = 1)
            #self.toolbar = NavigationToolbar2TkAgg(self.canvas, master)
            #self.toolbar.pack_configure(side=tk.BOTTOM)
        else:
            self.fig = plt.figure()
            self.canvas = self.fig.canvas
        self._blit = blit and self.fig.canvas.supports_blit
        self.canvas.mpl_connect('scroll_event', self.scale_view)
        self.ax_rect = self.fig.add_subplot(111, label = 'rect', xlim = (-5, 5), ylim=(-5, 5))
        self.ax_3d = self.fig.add_subplot(111, label = '3d', projection = '3d', zlim = (-1, 1))
        self.ax_polar = self.fig.add_subplot(111, label = 'polar', projection = 'polar')
        self.ax_main = self.ax_rect
        self.ax_x = self.fig.add_subplot(222, label = 'x', visible = False)
        self.ax_y = self.fig.add_subplot(223, label = 'y', visible = False, sharex = self.ax_x)
        self.ax_z = self.fig.add_subplot(224, label = 'z', visible = False, sharex = self.ax_x)
        self._3d = False
        self.mode = mode
        self.bg = dict()
        self.trajectories = []
        if master is None:
            self.fig.show()
        self.set_mode(mode)
        self.set_proj('rect')

    def get_diagonal(self):
        lim = self.get_limits()
        p1, p2 = np.array(lim).transpose()
        return helpers.distance(p1, p2)

    def get_limits(self):
        xlim = self.ax_main.get_xlim()
        ylim = self.ax_main.get_ylim()
        if not hasattr(self.ax_main, 'get_zlim'):
            return xlim, ylim, (-1.0, 1.0)
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

    def set_proj(self, projection = 'rect'):
        geom = self.ax_main.get_geometry()
        self._3d = False
        if projection == 'rect':
            self.ax_main = self.ax_rect
        elif projection == '3d':
            self.ax_main = self.ax_3d
            self._3d = True
        elif projection == 'polar':
            self.ax_main = self.ax_polar
        else:
            raise ValueError("Unknown projection: %s" % projection)
        for ax in self.ax_rect, self.ax_3d, self.ax_polar:
            if ax in self.fig.get_axes():
                ax.set_visible(False)
                self.fig.delaxes(ax)
        self.ax_main.set_visible(True)
        self.ax_main.change_geometry(*geom)
        self.fig.add_subplot(self.ax_main)
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
        self.mode = mode

    def clear(self, tmax = 1.0):
        xlim = self.ax_rect.get_xlim()
        ylim = self.ax_rect.get_ylim()
        zlim = self.ax_3d.get_zlim()
        tlim = 0.0, tmax
        for a in (self.ax_rect, self.ax_3d, self.ax_polar,
                self.ax_x, self.ax_y, self.ax_z):
            a.clear()
            a.grid(True)
        # The t-limit (time axis) is shared between ax_x, ax_y and ax_z
        self.ax_x.set_xlim(tlim)
        #self.ax_polar.set_xlim(self.ax_polar.get_xlim())
        self.ax_polar.set_ylim(self.ax_polar.get_ylim())
        self.set_limits(xlim, ylim, zlim)

    def set_limits(self, xlim, ylim, zlim):
        self.ax_rect.set_xlim(xlim)
        self.ax_rect.set_ylim(ylim)
        self.ax_3d.set_xlim(xlim)
        self.ax_3d.set_ylim(ylim)
        self.ax_3d.set_zlim(zlim)
        self.ax_x.set_ylim(xlim)
        self.ax_y.set_ylim(ylim)
        self.ax_z.set_ylim(zlim)

    def scale_view(self, evt):
        if evt.inaxes is self.ax_3d:
            return
        if not evt.inaxes is self.ax_main:
            return
        if evt.button == "up":
            scale = 0.8
        else:
            scale = 1.25
        pos = evt.xdata, evt.ydata
        xlim, ylim, zlim = self.get_limits()
        if evt.inaxes is self.ax_polar:
            rlim = 0.0, scale*ylim[1]
            self.ax_polar.set_ylim(rlim)
            return
        if evt.key is None or evt.key == 'x':
            xlim = helpers.scale_domain(xlim, scale, evt.xdata)
        if evt.key is None or evt.key == 'y':
            ylim = helpers.scale_domain(ylim, scale, evt.ydata)
        self.set_limits(xlim, ylim, zlim)
        self.draw(force = True)

    def draw(self, force = False):
        if force:
            self.fig.canvas.draw()
        if self._blit:
            self.fig.canvas.blit()
        else:
            self.fig.canvas.draw_idle()

    # NOTE: the draw_* functions must not call draw or blit. 
    def draw_nullclines(self, field, **kwargs):
        if self.ax_main is self.ax_3d:
            return
        x, y, u, v = field
        self.ax_main.contour(x, y, u, [0.0], **kwargs)
        self.ax_main.contour(x, y, v, [0.0], **kwargs)

    def draw_quiver(self, field, scaled = True, **kwargs):
        if self.ax_main is self.ax_3d:
            return
        x, y, u, v = field
        if scaled:
            mag = np.ma.masked_less(np.sqrt(u*u+v*v), 1.0e-4)
            u = u/mag
            v = v/mag
        self.ax_main.quiver(x, y, u, v, **kwargs)

    def draw_trajectory(self, traj, t_anim = None):
        if t_anim is None:
            for m in range(4):
                traj.marker[m].set_visible(False)
            traj.marker["start"].set_visible(True)
            traj.marker["3d"].set_visible(False)
            traj.arrow.set_visible(True)
            traj.line[0].set_xdata(traj.x)
            traj.line[0].set_ydata(traj.y)
            for l in 1, 2, 3:
                traj.line[l].set_xdata(traj.t)
                traj.line[l].set_ydata(traj.points[:,l-1])
            traj.line["3d"]._verts3d = (traj.x, traj.y, traj.z)
        else: # During animation
            if t_anim > traj.t[-1]:
                t_anim = traj.t[-1]
            last = traj.prev_index(t_anim)+1
            x, y, z = traj.x_ip(t_anim), traj.y_ip(t_anim), traj.z_ip(t_anim)
            tdata = np.append(traj.t[:last], t_anim)
            xdata = np.append(traj.x[:last], x)
            ydata = np.append(traj.y[:last], y)
            zdata = np.append(traj.z[:last], z)
            traj.line[0].set_xdata(xdata)
            traj.line[0].set_ydata(ydata)
            for l in 1, 2, 3:
                traj.line[l].set_xdata(tdata)
            traj.line[1].set_ydata(xdata)
            traj.line[2].set_ydata(ydata)
            traj.line[3].set_ydata(zdata)
            traj.marker[0].set_xdata([x])
            traj.marker[0].set_ydata([y])
            for m in 1, 2, 3:
                traj.marker[m].set_xdata([t_anim])
            traj.marker[1].set_ydata([x])
            traj.marker[2].set_ydata([y])
            traj.marker[3].set_ydata([z])
            traj.marker["3d"]._verts3d = ([x], [y], [z])
            for m in range(4):
                traj.marker[m].set_visible(True)
            traj.marker["start"].set_visible(True)
            traj.marker["3d"].set_visible(True)
        # Draw only the relevant artists
        if not self._3d:
            self._draw_artist(traj.line[0])
            self._draw_artist(traj.arrow)
            self._draw_artist(traj.marker["start"])
            self._draw_artist(traj.marker[0])
        else:
            self._draw_artist(traj.line["3d"])
            self._draw_artist(traj.marker["3d"])
        if self.mode == 4:
            for axis, axes in zip((1,2,3), (self.ax_x, self.ax_y, self.ax_z)):
                self._draw_artist(traj.line[axis], axes)
                self._draw_artist(traj.marker[axis], axes)

    def add_trajectory(self, traj, *args, **kwargs):
        self.trajectories.append(traj)
        self.init_trajectory(traj, *args, **kwargs)

    def init_trajectory(self, traj, style = 'r-', mfc = 'none', marker = 'bo'):
        """Creates the line, start marker, arrow and animation marker"""
        traj.line = dict()
        traj.marker = dict()
        traj.line[0], = self.ax_main.plot(traj.x, traj.y, style)
        traj.line[1], = self.ax_x.plot(traj.t, traj.x, style)
        traj.line[2], = self.ax_y.plot(traj.t, traj.y, style)
        traj.line[3], = self.ax_z.plot(traj.t, traj.z, style)
        traj.line["3d"], = self.ax_3d.plot(traj.x, traj.y, traj.z, style)
        traj.marker["start"], = self.ax_main.plot([traj.start[0]], [traj.start[1]], marker, mfc = mfc)
        traj.marker["3d"], = self.ax_3d.plot([traj.start[0]], [traj.start[1]], [traj.start[2]], marker)
        traj.marker[0], = self.ax_main.plot([traj.start[0]], [traj.start[1]], marker)
        traj.marker[1], = self.ax_x.plot([0.0], [traj.start[0]], marker)
        traj.marker[2], = self.ax_y.plot([0.0], [traj.start[1]], marker)
        traj.marker[3], = self.ax_z.plot([0.0], [traj.start[2]], marker)
        for m in range(4):
            traj.marker[m].set_visible(False)
        #traj.marker_anim.set_visible(False)
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
        self.restore()#self.ax_main)
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
