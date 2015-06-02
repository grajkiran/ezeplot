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
from widgets import AboutDialog

parser = OptionParser()
parser.set_defaults(loglevel = logging.WARNING)
parser.add_option("-b", "--blit", action = "store_true", default = True)
parser.add_option("-n", "--no-blit", action = "store_false", dest = "blit")
parser.add_option("-v", "--verbose", action = "store_const", dest = "loglevel",
        const = logging.INFO)
parser.add_option("-d", "--debug", action = "store_const", dest = "loglevel",
        const = logging.DEBUG)
opts, args = parser.parse_args()

logging.basicConfig(level = opts.loglevel)

system = DynamicSystem('y', '-x + mu * (1 - x*x)*y', params = dict(mu = 1.0))
logging.debug("%g: Creating root window..." % uptime.uptime())
root = tk.Tk()
root.title("Ezeplot - Dynamical systems visualization")
root.protocol('WM_DELETE_WINDOW', root.quit)
try:
    icon = tk.PhotoImage(file = os.path.join(os.path.dirname(__file__), 'icon-win.ppm'))
    root.tk.call("wm", "iconphoto", root._w, icon)
except:
    icon = None
logging.debug("%g: Creating application..." % uptime.uptime())
app = AppWindow(root, system, blit = opts.blit, icon = icon)
logging.debug("%g: Entering mainloop..." % uptime.uptime())
root.call('wm', 'attributes', '.', '-topmost', '1')
root.call('wm', 'attributes', '.', '-topmost', '0')
root.mainloop()
root.destroy()
