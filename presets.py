#!/usr/bin/python

systems = {
        "Linear System": {
            "x": "a*x + b*y",
            "y": "c*x + d*y",
            "z": "0",
            "params": {"a": 2, "b": 2, "c": -2, "d": -3}
            },
        "Simple Pendulum": {
            "x": "y",
            "y": "sin(x)",
            "z": "0"
            },
        "Linear Oscillator":   {
            "x": "y",
            "y": "-omega * x - c * y",
            "z": "0",
            "params": {"c": 0.5, "omega": 1.0}
            },
        "Van Der Pol Oscillator":  {
            "x": "y",
            "y": "-x + mu * (1 - x*x)*y",
            "z": "0",
            "params": {"mu": 1.0}
            },
        "Modified Van Der Pol":  {
            "x": "y - mu * (x**3/3 - x)",
            "y": "-x",
            "z": "0",
            "params": {"mu": 10.0}
            },
        "Glycolysis limit cycle":  {
            "x": "-x + a*y+x**2*y",
            "y": "b - a*y -x**2*y",
            "z": "0",
            "params": {"a": 0.05, "b": 0.5}
            },
        "Non isolated FP": {
            "x": "y",
            "y": "-2*mu*y - omega**2 * x",
            "z": "0",
            "params": {"mu": 1.0, "omega": 0.01}
            },
        "Glider problem": {
            "x": "-cos(x)/y + y",
            "y": "-sin(x) - D*y*y",
            "z": "0",
            "params": {"D": 0.0}
            },
        "Non linear center": {
            "x": "-y + a*x*(x*x + y*y)",
            "y": "x + a*y*(x*x + y*y)",
            "z": "0",
            "params": {"a": 0.0}
            },
        "Lorentz attractor": {
                "x": "sigma * (y-x)",
                "y": "x * (rho - z) - y",
                "z": "x*y - beta * z",
                "params": {"sigma": 10.0, "beta": 8.0/3.0, "rho": 28.0},
                "xlim": [-20.0, 20.0],
                "ylim": [-30.0, 30.0],
                "zlim": [0.0, 50],
                "tmax": 100.0,
                "projection": "3D"
                },
        "Duffing oscillator": {
                "x": "y",
                "y": "-d*y + x - x**3 + gamma*cos(w*z)",
                "z": "1.0",
                "xlim": [-10.0, 10.0],
                "ylim": [-10.0, 10.0],
                "zlim": [0, 50]
                }
        }
