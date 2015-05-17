#!/usr/bin/env python3
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
import uptime
try:
    import tkinter as tk
except ImportError:
    import Tkinter as tk
import matplotlib
matplotlib.use('TkAgg')
import sys
from optparse import OptionParser
import os.path
import logging

# Local imports
from plotting import Figure
from gui import AppWindow
from dynsystem import DynamicSystem


parser = OptionParser()
#parser.add_option("-e", "--embed", action = "store_true", default = True)
#parser.add_option("-w", "--windowed", action = "store_false", dest = 'embed')
parser.add_option("-b", "--blit", action = "store_true", default = True)
parser.add_option("-n", "--no-blit", action = "store_false", dest = "blit")
opts, args = parser.parse_args()

system = DynamicSystem('y', '-x + mu * (1 - x*x)*y', params = dict(mu = 1.0))
#fig = Figure(blit = opts.blit)
logging.debug(uptime.uptime(), "Creating root window...")
root = tk.Tk()
root.title("EzePlot - Dynamical systems visualization")
root.protocol('WM_DELETE_WINDOW', root.quit)
root.tk.call("wm", "iconphoto", root._w, tk.PhotoImage(file = 'icon.png'))
#root.attributes('-fullscreen', True)
#root.attributes('-zoomed', True)
#fig.bind('close_event', lambda evt: root.quit())
logging.debug(uptime.uptime(), "Creating application...")
app = AppWindow(root, system, blit = opts.blit)#, embedded = opts.embed)
# TODO: Set the geometry of the control window and figure windows
logging.debug(uptime.uptime(), "Entering mainloop...")
root.mainloop()
