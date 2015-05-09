#!/usr/bin/python3
# External libraries
try:
    import tkinter as tk
except ImportError:
    import Tkinter as tk
import matplotlib
matplotlib.use('TkAgg')
import sys
from optparse import OptionParser

# Local imports
from plotting import Figure
from gui import AppWindow
from dynsystem import DynamicSystem
import uptime


parser = OptionParser()
#parser.add_option("-e", "--embed", action = "store_true", default = True)
#parser.add_option("-w", "--windowed", action = "store_false", dest = 'embed')
parser.add_option("-b", "--blit", action = "store_true", default = True)
parser.add_option("-n", "--no-blit", action = "store_false", dest = "blit")
opts, args = parser.parse_args()

system = DynamicSystem('y', '-x + mu * (1 - x*x)*y', params = dict(mu = 1.0))
#fig = Figure(blit = opts.blit)
print(uptime.uptime(), "Creating root window...")
root = tk.Tk()
root.title("EzePlot - Dynamic systems visualization")
root.protocol('WM_DELETE_WINDOW', root.quit)
#root.attributes('-fullscreen', True)
#root.attributes('-zoomed', True)
#fig.bind('close_event', lambda evt: root.quit())
print(uptime.uptime(), "Creating application...")
app = AppWindow(root, system, blit = opts.blit)#, embedded = opts.embed)
# TODO: Set the geometry of the control window and figure windows
print(uptime.uptime(), "Entering mainloop...")
root.mainloop()
