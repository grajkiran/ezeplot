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
except:
    import Tkinter as tk
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import logging
import numpy as np
import helpers

import uptime
matplotlib.rcParams['toolbar'] = 'None'
matplotlib.rc('font', family = 'serif')
#matplotlib.rc('text', usetex = True)
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
#matplotlib.rc('font', size = 16, weight = 'bold')

class Figure:
    def __init__(self, master, temporal = False, blit = True):
        self.fig = matplotlib.figure.Figure()
        self.canvas = FigureCanvasTkAgg(self.fig, master = master)
        #self.canvas.get_tk_widget().pack(side = tk.LEFT, fill = tk.BOTH, expand = 1)
        self._blit = blit and self.fig.canvas.supports_blit
        #self.canvas.mpl_connect('scroll_event', self.scale_view)
        logging.debug("%g: Creating subplots main..." % uptime.uptime())
        self.ax_rect = self.fig.add_subplot(111, label = 'rect', xlim = (-5, 5), ylim=(-5, 5))
        self.ax_3d = self.fig.add_subplot(111, label = '3d', projection = '3d', zlim = (-1, 1))
        self.ax_polar = self.fig.add_subplot(111, label = 'polar', projection = 'polar')
        self.ax_main = self.ax_rect
        logging.debug("%g: Creating subplots time series..." % uptime.uptime())
        self.ax_x = self.fig.add_subplot(222, label = 'x', visible = False)
        self.ax_y = self.fig.add_subplot(223, label = 'y', visible = False, sharex = self.ax_x)
        self.ax_z = self.fig.add_subplot(224, label = 'z', visible = False, sharex = self.ax_x)
        self._3d = False
        self.temporal = temporal
        self.bg = dict()
        #self.trajectories = []
        logging.debug("%g: Setting mode..." % uptime.uptime())
        #self.set_mode(temporal)
        logging.debug("%g: Setting projection..." % uptime.uptime())
        self.set_proj('rect', temporal)

    def get_diagonal(self):
        lim = self.get_limits()
        p1, p2 = np.array(lim).transpose()
        return helpers.distance(p1, p2)

    def get_limits(self):
        xlim = self.ax_main.get_xlim()
        ylim = self.ax_main.get_ylim()
        if not hasattr(self.ax_main, 'get_zlim'):
            return xlim, ylim, (-5.0, 5.0)
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
        if not ax in self.fig.get_axes():
            return
        if not artist in ax.get_children():
            ax.add_artist(artist)
        if self._blit:
            ax.draw_artist(artist)

    def set_proj(self, projection, temporal):
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
        self.fig.add_axes(self.ax_main)
        self.set_mode(temporal)
        #self.ax_main.change_geometry(*geom)
        self.draw(force = True)

    def set_mode(self, temporal):
        for ax in self.ax_x, self.ax_y, self.ax_z:
            ax.set_visible(False)
            if ax in self.fig.get_axes():
                self.fig.delaxes(ax)
        if not temporal:
            self.ax_main.change_geometry(1, 1, 1)
        else:
            self.ax_main.change_geometry(2, 2, 1)
            for ax in self.ax_x, self.ax_y:
                self.fig.add_axes(ax)
                ax.set_visible(True)
            if self._3d:
                self.fig.add_axes(self.ax_z)
                self.ax_z.set_visible(True)
        self.temporal = temporal

    def __set_3d_mode(self):
        self.ax_3d.mouse_init(zoom_btn = [], rotate_btn = [1])
        self.ax_3d.grid(False)
        self.ax_3d.set_frame_on(False)
        #self.ax_3d.set_axis_off()
        #self.ax_3d.set_xticks([])
        #self.ax_3d.set_yticks([])
        #self.ax_3d.set_zticks([])
    def clear(self, tmax = 1.0):
        xlim = self.ax_rect.get_xlim()
        ylim = self.ax_rect.get_ylim()
        zlim = self.ax_3d.get_zlim()
        var_x = 'X'
        var_y = 'Y'
        if self.ax_main is self.ax_polar:
            var_x = r'\theta'
            var_y = 'r'
        title_x = "Time series $(%s)$" % var_x
        title_y = "Time series $(%s)$" % var_y
        title_z = "Time series $(Z)$"
        tlim = 0.0, tmax
        for a in (self.ax_rect, self.ax_3d, self.ax_polar,
                self.ax_x, self.ax_y, self.ax_z):
            a.clear()
            a.grid(True)
        self.fig.subplots_adjust(top = 0.85, left = 0.1, right = 0.9,
                bottom = 0.1, hspace = 0.4)
        self.__set_3d_mode()
        self.ax_rect.text(0.5, 1.1, "Phase portrait $(XY)$", size = 16, weight = 'bold',
                transform = self.ax_main.transAxes, ha = 'center')
        self.ax_polar.text(0.5, 1.1, r"Phase portrait $(\theta r)$", size = 16, weight = 'bold',
                transform = self.ax_polar.transAxes, ha = 'center')
        self.ax_3d.text2D(0.5, 1.1, "Phase portrait $(XYZ)$", size = 16, weight = 'bold',
                transform = self.ax_3d.transAxes, ha = 'center')
        for ax, title in zip((self.ax_x, self.ax_y, self.ax_z), (title_x, title_y, title_z)):
            ax.text(0.5, 1.1, title, size = 16, weight = 'bold',
                transform = ax.transAxes, ha = 'center')
        #self.ax_x.text(0.5, 1.1, "Time series $(X)$", size = 16, weight = 'bold',
        #        transform = self.ax_x.transAxes, ha = 'center')
        #self.ax_y.text(0.5, 1.1, "Time series $(Y)$", size = 16, weight = 'bold',
        #        transform = self.ax_y.transAxes, ha = 'center')
        #self.ax_z.text(0.5, 1.1, "Time series $(Z)$", size = 16, weight = 'bold',
        #        transform = self.ax_z.transAxes, ha = 'center')
        #polar_message = "$x=\\theta$\n$y=r$"
        #self.ax_polar.annotate(polar_message, (-40, -75), size = 16, color = 'blue',
        #        xycoords = 'figure points')
        self.ax_rect.set_xlabel(r'X')
        self.ax_rect.set_ylabel(r'Y')
        self.ax_x.set_xlabel(r'Time')
        self.ax_y.set_xlabel(r'Time')
        self.ax_z.set_xlabel(r'Time')
        self.ax_x.set_ylabel(r'$%s$' % var_x)
        self.ax_y.set_ylabel(r'$%s$' % var_y)
        self.ax_z.set_ylabel(r'Z')
        self.ax_3d.set_xlabel(r'X')
        self.ax_3d.set_ylabel(r'Y')
        self.ax_3d.set_zlabel(r'Z')
        # The t-limit (time axis) is shared between ax_x, ax_y and ax_z
        self.ax_x.set_xlim(tlim)
        #self.ax_polar.set_xlim(self.ax_polar.get_xlim())
        self.ax_polar.set_ylim(self.ax_polar.get_ylim())
        self.set_limits(xlim, ylim, zlim)

    def set_limits(self, xlim, ylim, zlim):
        self.ax_rect.set_xlim(xlim)
        self.ax_rect.set_ylim(ylim)
        if max(ylim) > 0:
            self.ax_polar.set_rlim([0, max(ylim)])
        self.ax_3d.set_xlim(xlim)
        self.ax_3d.set_ylim(ylim)
        self.ax_3d.set_zlim(zlim)
        self.ax_x.set_ylim(xlim)
        self.ax_y.set_ylim(ylim)
        self.ax_z.set_ylim(zlim)
        if self.ax_main is self.ax_polar:
            self.ax_x.set_ylim((0, 2*np.pi))
            self.ax_y.set_ylim(self.ax_polar.get_ylim())

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
        else:
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

    def draw_quiver3d(self, field, scaled = True, **kwargs):
        if self.ax_main is not self.ax_3d:
            return
        x, y, z, u, v, w= field
        if scaled:
            mag = np.ma.masked_less(np.sqrt(u*u+v*v+w*w), 1.0e-4)
            u = u/mag
            v = v/mag
            w = w/mag
        self.ax_main.quiver(x, y, z, u, v, w)

    def draw_quiver(self, field, scaled = True, **kwargs):
        if self.ax_main is self.ax_3d:
            return
        x, y, u, v = field
        if scaled:
            mag = np.ma.masked_less(np.sqrt(u*u+v*v), 1.0e-4)
            u = u/mag
            v = v/mag
        self.ax_main.quiver(x, y, u, v, **kwargs)

    def draw_fp(self, *fixed_points):
        if len(fixed_points) == 0:
            return
        fps = np.array(fixed_points)
        if self.ax_main is self.ax_3d:
            art, = self.ax_main.plot(fps[:,0], fps[:,1], fps[:,2], 'r*', markersize = 8, mec = 'red')
        else:
            art, = self.ax_main.plot(fps[:,0], fps[:,1], 'r*', markersize = 8, mec='red')
        self._draw_artist(art)

    def _draw_trajectory(self, traj, t_anim = None):
        #logging.debug("Drawing trajectory: %s" % str(traj))
        if t_anim is None:
            for m in range(4):
                traj.marker[m].set_visible(False)
            traj.line[0].set_xdata(traj.x)
            traj.line[0].set_ydata(traj.y)
            if self._3d:
                traj.line[0].set_3d_properties(traj.z)
            for l in 1, 2, 3:
                traj.line[l].set_xdata(traj.t)
                traj.line[l].set_ydata(traj.points[:,l-1])
            if self.ax_main is self.ax_polar:
                traj.line[1].set_ydata(traj.points[:,0]%(2*np.pi))
        else: # During animation
            if t_anim > traj.t[-1]:
                t_anim = traj.t[-1]
            last = traj.prev_index(t_anim)#+1
            x, y, z = traj.points[last]
            tdata = np.append(traj.t[:last], t_anim)
            xdata = np.append(traj.x[:last], x)
            ydata = np.append(traj.y[:last], y)
            zdata = np.append(traj.z[:last], z)
            traj.line[0].set_xdata(xdata)
            traj.line[0].set_ydata(ydata)
            if self._3d:
                traj.line[0].set_3d_properties(zdata)
            for l in 1, 2, 3:
                traj.line[l].set_xdata(tdata)
            if self.ax_main is self.ax_polar:
                traj.line[1].set_ydata(xdata%(2*np.pi))
            else:
                traj.line[1].set_ydata(xdata)
            traj.line[2].set_ydata(ydata)
            traj.line[3].set_ydata(zdata)
            traj.marker[0].set_xdata([x])
            traj.marker[0].set_ydata([y])
            for m in 1, 2, 3:
                traj.marker[m].set_xdata([t_anim])
            if self.ax_main is self.ax_polar:
                traj.marker[1].set_ydata([x%(2*np.pi)])
            else:
                traj.marker[1].set_ydata([x])
            traj.marker[2].set_ydata([y])
            traj.marker[3].set_ydata([z])
            if self._3d:
                traj.marker[0].set_3d_properties([z])
            for m in range(4):
                traj.marker[m].set_visible(True)
        # Draw only the relevant artists
        self._draw_artist(traj.line[0])
        self._draw_artist(traj.marker[0])
        if self.temporal:
            for axis, axes in zip((1,2,3), (self.ax_x, self.ax_y, self.ax_z)):
                self._draw_artist(traj.line[axis], axes)
                self._draw_artist(traj.marker[axis], axes)

    #def add_trajectory(self, traj, *args, **kwargs):
    #    #self.trajectories.append(traj)
    #    self.init_trajectory(traj, *args, **kwargs)

    def add_trajectory(self, traj, style = 'b-', mfc = 'none', marker = 'ro'):
        """Creates the line, start marker, arrow and animation marker"""
        if hasattr(traj, 'style'):
            style = traj.style
        traj.line = dict()
        traj.marker = dict()
        if self._3d:
            traj.line[0], = self.ax_main.plot(traj.x, traj.y, traj.z, style)
            traj.marker[0], = self.ax_main.plot([traj.start[0]], [traj.start[1]], [traj.start[2]], marker)
        else:
            traj.line[0], = self.ax_main.plot(traj.x, traj.y, style)
            traj.marker[0], = self.ax_main.plot([traj.start[0]], [traj.start[1]], marker)
        traj.line[1], = self.ax_x.plot(traj.t, traj.x, style)
        traj.line[2], = self.ax_y.plot(traj.t, traj.y, style)
        traj.line[3], = self.ax_z.plot(traj.t, traj.z, style)
        traj.marker[1], = self.ax_x.plot([0.0], [traj.start[0]], marker)
        traj.marker[2], = self.ax_y.plot([0.0], [traj.start[1]], marker)
        traj.marker[3], = self.ax_z.plot([0.0], [traj.start[2]], marker)
        for m in range(4):
            traj.marker[m].set_visible(False)
        self._draw_trajectory(traj)

    def anim_update(self, time, trajectories):
        self.restore()#self.ax_main)
        for traj in trajectories:
            self._draw_trajectory(traj, time)
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
