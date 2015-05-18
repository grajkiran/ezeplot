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
from collections import OrderedDict

systems = OrderedDict()
systems["Simple pendulum"] = {
            "x": "y",
            "y": "sin(x)",
            "z": "0",
            "locations": [(0.5, 0.5), (0.5, -0.5), (2, 1), (-2, 1)],
            "xlim": [-5, 5],
            "ylim": [-5, 5],
            "tmax": 25.0,
            "projection": "2D",
            }
systems["Linear damped oscillator"] = {
            "x": "y",
            "y": "-omega * x - c * y",
            "z": "0",
            "params": {"c": 0.5, "omega": 1.0},
            "locations": [(4, 2), (-4, -2)],
            "xlim": [-5, 5],
            "ylim": [-5, 5],
            "tmax": 25.0,
            "projection": "2D",
            }
systems["Van der Pol oscillator"] = {
            "x": "y",
            "y": "-x + mu * (1 - x*x)*y",
            "z": "0",
            "params": {"mu": 2.5},
            "locations": [(0.1, 0.1), (4.0, 4.0), (-4, -4), (4, -4), (-4, 4)],
            "xlim": [-5, 5],
            "ylim": [-5, 5],
            "tmax": 25.0,
            "projection": "2D",
            }
systems["Belousov-Zhabotinski"] = {
            "x": "a - x - 4*x*y/(1+x*x)",
            "y": "b*x *(1- y/(1+x*x))",
            "z": "0",
            "params": {"a": 10, "b": 1},
            "xlim": [0, 5],
            "ylim": [0, 10],
            "tmax": 50,
            "projection": "2D",
            "locations": [(4,4), (2,6)],
            }
systems["Lorentz attractor"] = {
            "x": "sigma * (y-x)",
            "y": "x * (rho - z) - y",
            "z": "x*y - beta * z",
            "params": {"sigma": 10.0, "beta": 8.0/3.0, "rho": 28.0},
            "xlim": [-20.0, 20.0],
            "ylim": [-30.0, 30.0],
            "zlim": [0.0, 50],
            "tmax": 100.0,
            "projection": "3D",
            "locations": [(0.5, 0.5)],
            "reverse": False,
            }
systems["Rossler attractor"] = {
            "x": "-y - z",
            "y": "x + a*y",
            "z": "b + z*(x - c)",
            "params": {"a": 0.1, "b": 0.1, "c": 14},
            "xlim": [-25.0, 25.0],
            "ylim": [-25.0, 25.0],
            "zlim": [-10.0, 40.0],
            "tmax": 100.0,
            "projection": "3D",
            "locations": [(-10.5, 10.5), (-10.5, 10.6)],
            }
systems["User defined"] = {
            "x": "0",
            "y": "0",
            "z": "0",
            "xlim": [-5, 5],
            "ylim": [-5, 5],
            "tmax": 25.0,
            "dt": 0.05,
            "projection": "2D",
            }
#systems["Duffing oscillator"] = {
#                "x": "y",
#                "y": "-d*y + x - x**3 + gamma*cos(w*z)",
#                "z": "1.0",
#                "xlim": [-5.0, 5.0],
#                "ylim": [-5.0, 5.0],
#                "zlim": [0, 50],
#                "params": {"d": 0, "gamma": -2, "w": 1},
#                "locations": [(0.7, -0.3, 0), (-0.18, 4.36, 0)],
#                "reverse": False,
#                "tmax": 50,
#                }
