#!/usr/bin/python3
import numpy as np
from scipy.integrate import ode
from scipy.interpolate import interp1d
from math import isinf
import helpers
from numpy import *
abs = absolute

MAX_CHARS = 6

class Trajectory:
    def __init__(self, data, startidx = 0):
        self.startidx = startidx
        self.set_data(data, startidx)

    def set_data(self, data, startidx = None):
        if startidx is None:
            startidx = self.startidx
        self.startidx = startidx
        self.data = np.array(data)
        self.time = self.data[:,0]
        self.dist = self.data[:,1]
        self.points = self.data[:,2:]
        self.t = self.time
        self.x = self.points[:,0]
        self.y = self.points[:,1]
        self.z = self.points[:,2]
        self.start = self.points[startidx]
        self.length = self.dist[-1] - self.dist[0]

    def prev_index(self, t):
        """Returns the maximum index of the datapoint whose timestamp is not
        more than t"""
        if t < self.t[0] or t > self.t[-1]:
            raise ValueError("t is out of bounds")
        prev = 0
        next = len(self.t)-1
        if t == self.t[prev]:
            return prev
        if t == self.t[next]:
            return next
        while(next - prev > 1):
            mid = (next+prev)//2
            if self.t[mid] == t:
                return mid
            if t > self.t[mid]:
                prev = mid
            else:
                next = mid
        return prev

class DynamicSystem:
    """A 3D dynamic system"""
    def __init__(self, x_dot = "x", y_dot = "y", z_dot = "0.0", params = dict()):
        self.__x_dot = str(x_dot)
        self.__y_dot = str(y_dot)
        self.__z_dot = str(z_dot)
        self.params = dict(params)
        self.limits = [None, None, None]
        self.periodic = [False, False, False]
        self.update_handlers = []
        self.update()

    def add_update_handler(self, handler):
        if callable(handler):
            self.update_handlers.append(handler)

    def get(self, which = None):
        if which == 1 or which == 'x':
            return self.__x_dot
        elif which == 2 or which == 'y':
            return self.__y_dot
        elif which == 3 or which == 'z':
            return self.__z_dot
        return self.__x_dot, self.__y_dot, self.__z_dot

    def update(self, x_dot = None, y_dot = None, z_dot = None):
        if x_dot is None:
            x_dot = self.__x_dot
        if y_dot is None:
            y_dot = self.__y_dot
        if z_dot is None:
            z_dot = self.__z_dot
        x_dot = self.__sanitize(x_dot)
        y_dot = self.__sanitize(y_dot)
        z_dot = self.__sanitize(z_dot)
        self.__code_x = compile(x_dot, '<string>', 'eval')
        self.__code_y = compile(y_dot, '<string>', 'eval')
        self.__code_z = compile(z_dot, '<string>', 'eval')
        self.__x_dot = x_dot
        self.__y_dot = y_dot
        self.__z_dot = z_dot
        self.__params_update()
        for p in self.params:
            if len(p) > MAX_CHARS:
                raise ValueError("Parameter name exceeds %d characters" % MAX_CHARS)
        self.__call__([0.8934,0.6572,0.8723])
        for handler in self.update_handlers:
            handler()

    def __sanitize(self, text):
        if '"' in (text) or "'" in (text):
            raise ValueError("Invalid characters")
        return text.replace("^", "**")

    def __params_update(self):
        #from numpy import sqrt, sin, cos, tan, abs, exp, pi
        x, y, z = 0.0, 0.0, 0.0
        req_params = set(self.__code_x.co_names)
        req_params.update(set(self.__code_y.co_names))
        req_params.update(set(self.__code_z.co_names))
        req_params.difference_update(globals().keys())
        req_params.difference_update(locals().keys())
        for p in req_params:
            if not p in self.params:
                self.params[p] = 1.0
        for p in list(self.params.keys()):
            if not p in req_params:
                del self.params[p]

    @staticmethod
    def normalize_coord(c, period):
        if period is None:
            return c
        low, high = sorted(period)
        mag = high - low
        c = np.array(c, dtype='float64')
        for x in np.nditer(c, op_flags=['readwrite']):
            while x < low:
                x += mag
            while x > high:
                x -= mag
        return c

    def __call__(self, x, params = None):
        #from numpy import sqrt, sin, cos, tan, abs, exp, pi, log, arctan
        f_params = globals()
        size = len(x)
        if size == 2:
            x, y = np.array(x)
            z = np.zeros_like(x)
        else:
            x, y, z = np.array(x)
        if self.periodic[0]:
            x = self.normalize_coord(x, self.limits[0])
        if self.periodic[1]:
            y = self.normalize_coord(y, self.limits[1])
        if self.periodic[2]:
            z = self.normalize_coord(z, self.limits[2])
        # We first initialize x_dot, y_dot, z_dot to 1.0 and then scale them with the
        # evaulated expression. This ensures that the resulting *_dot 
        # are always numpy arrays and don't get clobbered if the user supplied
        # expression evaulates to a constant.
        x_dot = np.ones_like(x)
        y_dot = np.ones_like(y)
        z_dot = np.ones_like(z)
        f_params.update(locals())
        x_dot *= eval(self.__code_x, f_params, params or self.params)
        y_dot *= eval(self.__code_y, f_params, params or self.params)
        z_dot *= eval(self.__code_z, f_params, params or self.params)
        return [x_dot, y_dot, z_dot]

    def trajectory(self, start, tmax, limits = None, threshold = 0.0,
            bidirectional = True, **kwargs):
        t_forw = self.__compute_trajectory(start, tmax,
                threshold, **kwargs)
        if not bidirectional:
            return Trajectory(t_forw, 0)
        # Remove the first element of the reversed trajectory as it is the
        # start point that is already present in the forward trajectory.
        t_back = self.__compute_trajectory(start, -tmax,
                threshold, **kwargs)[1:]
        t_back.reverse()
        return Trajectory(t_back+t_forw, len(t_back))

    def __compute_trajectory(self, start, tmax, threshold, **kwargs):
        traj = []
        threshold = abs(threshold)
        def accumulate(t, pos):
            pos = list(pos)
            pos_norm = list(pos)
            if len(traj) == 0:
                traj.append([t, 0.0] + pos_norm)
                return 0
            s_old = traj[-1][1]
            pos_old = traj[-1][2:]
            s = helpers.distance(pos_old, pos)
            if tmax > 0:
                s_cum = s_old + s
            else:
                s_cum = s_old - s
            traj.append([t, s_cum] + pos_norm)
            if not helpers.is_inside(pos_norm, self.limits):
                print("Crossed domain - stopping integration.")
                return -1
            # Check if we are within threshold of a fixed point.
            # We are if our velocity is decreasing and we have not moved much.
            v_old = helpers.mag(self(pos_old))
            v = helpers.mag(self(pos_norm))
            if s < threshold and v < v_old and len(traj) > 10:
                print("No significant change - stopping integration.")
                return -1
            return 0
        rk4 = ode(lambda t, pos: self(pos)).set_integrator('dopri5', **kwargs)
        rk4.set_solout(accumulate)
        rk4.set_initial_value(start, 0.0)
        rk4.integrate(tmax)
        #print("Integration took %g seconds" % t.seconds())
        return traj

def main():
    d = DynamicSystem()
    p = [-5, 5]

if __name__ == '__main__':
    main()
