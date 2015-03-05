#!/usr/bin/python3
import numpy as np
from scipy.integrate import ode
from scipy.interpolate import interp1d
from math import isinf
import helpers

def is_inside(p, limits):
    if limits is None:
        return True
    for i in range(len(p)):
        lower = min(limits[i])
        upper = max(limits[i])
        if p[i] < lower or p[i] > upper:
            return False
    return True

class Trajectory:
    def __init__(self, data, startidx = 0):
        self.startidx = startidx
        self.data = np.array(data)
        self.time = self.data[:,0]
        self.dist = self.data[:,1]
        self.points = self.data[:,2:]
        self.t = self.time
        self.x = self.points[:,0]
        self.y = self.points[:,1]
        self.s_ip = interp1d(self.time, self.dist)
        self.t_ip = interp1d(self.dist, self.time)
        self.x_ip = interp1d(self.time, self.x)
        self.y_ip = interp1d(self.time, self.y)
        self.start = self.points[startidx]
        self.length = self.dist[-1] - self.dist[0]

#class TrajAccumulator(list):
#    def __init__(self, limits = None, threshold = 0.0):
#        self.limits = limits
#        self.threshold = abs(threshold)
#
#    def accumulate(self, t, pos):
#        pos = list(pos)
#        if len(self) == 0:
#            self.append([t, 0.0] + pos)
#            return 0
#        pos_old = self[-1][2:]
#        s = helpers.distance(pos_old, pos)
#        self.append([t, s] + pos)
#        if s < self.threshold or not is_inside(pos, self.limits):
#            return -1
#        return 0

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

    def trajectory(self, start, tmax, limits = None, threshold = 0.0,
            bidirectional = True, **kwargs):
        t_forw = self.__compute_trajectory(start, tmax, limits,
                threshold, **kwargs)
        if not bidirectional:
            return Trajectory(t_forw, 0)
        # Remove the first element of the reversed trajectory as it is the
        # start point that is already present in the forward trajectory.
        t_back = self.__compute_trajectory(start, -tmax, limits,
                threshold, **kwargs)[1:]
        t_back.reverse()
        return Trajectory(t_back+t_forw, len(t_back))

    def __compute_trajectory(self, start, tmax, limits, threshold, **kwargs):
        traj = []
        threshold = abs(threshold)
        def accumulate(t, pos):
            pos = list(pos)
            if len(traj) == 0:
                traj.append([t, 0.0] + pos)
                return 0
            s_old = traj[-1][1]
            pos_old = traj[-1][2:]
            s = helpers.distance(pos_old, pos)
            if tmax > 0:
                s_cum = s_old + s
            else:
                s_cum = s_old - s
            traj.append([t, s_cum] + pos)
            if s < threshold or not is_inside(pos, limits):
                return -1
            return 0
        #tdata = TrajAccumulator(limits, threshold)
        rk4 = ode(lambda t, pos: self(pos)).set_integrator('dopri5', **kwargs)
        rk4.set_solout(accumulate)
        rk4.set_initial_value(start, 0.0)
        rk4.integrate(tmax)
        return traj

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
