#!/usr/bin/python3
import numpy as np
from scipy.integrate import ode
from math import isinf
import helpers

class Trajectory:
    def __init__(self, data, start = None):
        self.t = data[:,0]
        self.x = data[:,1]
        self.y = data[:,2]
        self.xy = data[:,1:]
        self.data = data
        self.start = start
        if start is None:
            self.start = self.data[0]
        self.direction = 0
        if self.start == self.data[0]:
            self.direction = 1
        elif self.start == self.data[-1]:
            self.direction = -1
        self.freq_x, self.amp_x = helpers.freq_domain(self.t, self.x)
        self.freq_y, self.amp_y = helpers.freq_domain(self.t, self.y)

class DynamicSystem:
    """A 2D dynamic system"""
    def __init__(self, x_dot = "x", y_dot = "y", params = dict()):
        self.__x_dot = x_dot
        self.__y_dot = y_dot
        self.params = dict(params)
        self.update_handlers = []
        self.update()

    def add_update_handler(self, handler):
        if callable(handler):
            self.update_handlers.append(handler)

    def get(self, which = None):
        if which == 1:
            return self.__x_dot
        elif which == 2:
            return self.__y_dot
        return self.__x_dot, self.__y_dot

    def update(self, x_dot = None, y_dot = None):
        if x_dot is None:
            x_dot = self.__x_dot
        if y_dot is None:
            y_dot = self.__y_dot
        if '"' in x_dot+y_dot or "'" in x_dot+y_dot: return False
        self.__code_x = compile(x_dot, '<string>', 'eval')
        self.__code_y = compile(y_dot, '<string>', 'eval')
        self.__x_dot = x_dot
        self.__y_dot = y_dot
        self.__params_update()
        self.__call__([0.8934,0.6572])
        for handler in self.update_handlers:
            handler()
            #handler(x_dot, y_dot, self.params)

    def __params_update(self):
        from numpy import sqrt, sin, cos, tan, abs, exp, pi
        x, y = 0.0, 0.0
        req_params = set(self.__code_x.co_names)
        req_params.update(set(self.__code_y.co_names))
        req_params.difference_update(locals().keys())
        for p in req_params:
            if not p in self.params:
                self.params[p] = 1.0
        for p in list(self.params.keys()):
            if not p in req_params:
                del self.params[p]
        return req_params

    def __call__(self, x, code_x = None, code_y = None, params = None):
        from numpy import sqrt, sin, cos, tan, abs, exp, pi
        x, y = x
        x_dot = eval(code_x or self.__code_x, locals(), params or self.params)
        y_dot = eval(code_y or self.__code_y, locals(), params or self.params)
        return [x_dot, y_dot]

    def trajectory(self, start, direction = 1, **kwargs):
        if 'reverse' in kwargs:
            if kwargs.pop('reverse'):
                direction = -1
        if direction > 0:
            t = self.__compute_trajectory(start, **kwargs)
        elif direction < 0:
            t = self.__compute_trajectory(start, reverse = True, **kwargs)
        else:
            t1 = self.__compute_trajectory(start, reverse = False, **kwargs)
            t2 = self.__compute_trajectory(start, reverse = True, **kwargs)
            t = np.append(t2, t1, 0)
        return Trajectory(t, start)

    def __compute_trajectory(self, start, dt = 0.01, max_steps = 100, explicit = False,
            xlim = None, ylim = None, ds_min = 1.0e-4, reverse = False, outer = 5):
        def is_inside(v, lim):
            if lim is None: return True
            a, b = sorted(lim)
            return v >= a and v <= b
        if reverse: dt *= -1.0
        traj = np.zeros((max_steps, 3))
        traj[0] = 0.0, start[0], start[1]
        rk4int = ode(lambda t, pos: self(pos)).set_integrator('dopri5')
        rk4int.set_initial_value(start, 0.0)
        ts_break = max_steps
        for tstep in range(1, max_steps):
            pos = np.array(traj[tstep-1][1:])
            time = tstep * dt
            if explicit:
                x, y = pos + dt * np.array(self(pos))
            else:
                x, y = rk4int.integrate(time)
                traj[tstep] = time, x, y
            # Break out if the trajectory crossed the limits more than 'outer'
            # steps ago.
            if is_inside(x, xlim) and is_inside(y, ylim):
                ts_break = max_steps # We are inside the domain
            else:
                if tstep - ts_break > outer:
                    break
                # Set ts_break to the current tstep only if not already set
                elif ts_break == max_steps:
                    ts_break = tstep
            # If the movement is very small, we must have reached a fixed point.
            if helpers.distance(pos, (x, y)) < ds_min: break
        if reverse:
            return traj[tstep::-1]
        return traj[:tstep+1]

presets = {
        'Linear System': (
            'a*x + b*y',
            'c*x + d*y',
            dict(a = 2, b = 2, c = -2, d = -3)),
        'Simple Pendulum': (
            'y',
            'sin(x)',
            dict()),
        'Linear Oscillator':   (
            'y',
            '-omega * x - c * y',
            dict(c = 0.5, omega = 1.0)),
        'Van Der Pol Oscillator':  (
            'y',
            '-x + mu * (1 - x*x)*y',
            dict(mu = 1.0)),
        'Modified Van Der Pol':  (
            'y - mu * (x**3/3 - x)',
            '-x',
            dict(mu = 10.0)),
        'Glycolysis limit cycle':  (
            '-x + a*y+x**2*y',
            'b - a*y -x**2*y',
            dict(a = 0.05, b = 0.5)),
        'Non isolated FP': (
            'y',
            '-2*mu*y - omega**2 * x',
            dict(mu = 1.0, omega = 0.01)),
        'Glider problem': (
            '-cos(x)/y + y',
            '-sin(x) - D*y*y',
            dict(D = 0.0)),
        'Non linear center': (
            '-y + a*x*(x*x + y*y)',
            'x + a*y*(x*x + y*y)',
            dict(a = 0.0)),
        }

def main():
    d = DynamicSystem(x_dot = "sin(2)*x+mu*y+0", params = dict(mu = 0.5))
    print(d.params.keys())
    l = d.trajectory([1.4, 1.5], xlim = (-5,5), ylim = (-5,5), reverse = True)
    print(l)
#    import IPython
#    IPython.embed()

if __name__ == '__main__':
    main()
