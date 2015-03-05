#!/usr/bin/python3
# External libraries
try:
    import tkinter as tk
except ImportError:
    import Tkinter as tk

# Local imports
from plotting import FigureWindow
from gui import AppWindow
from dynsystem import DynamicSystem

system = DynamicSystem('y', '-x + mu * (1 - x*x)*y', dict(mu = 1.0))
fig = FigureWindow()
root = tk.Tk()
root.protocol('WM_DELETE_WINDOW', root.quit)
fig.bind('close_event', lambda evt: root.quit())
app = AppWindow(root, system, fig)
# TODO: Set the geometry of the control window and figure windows
root.mainloop()
