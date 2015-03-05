#!/usr/bin/python3
# External libraries
try:
    import tkinter as tk
except ImportError:
    import Tkinter as tk
import sys
from optparse import OptionParser

# Local imports
from plotting import FigureWindow
from gui import AppWindow
from dynsystem import DynamicSystem


parser = OptionParser()
parser.add_option("-b", "--blit", action = "store_true", default = True)
parser.add_option("-n", "--no-blit", action = "store_false", dest = "blit")
opts, args = parser.parse_args()

system = DynamicSystem('y', '-x + mu * (1 - x*x)*y', dict(mu = 1.0))
fig = FigureWindow(blit = opts.blit)
root = tk.Tk()
root.protocol('WM_DELETE_WINDOW', root.quit)
fig.bind('close_event', lambda evt: root.quit())
app = AppWindow(root, system, fig)
# TODO: Set the geometry of the control window and figure windows
root.mainloop()
