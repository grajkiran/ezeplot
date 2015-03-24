#!/usr/bin/python

def sign(x):
    if x > 0: return 1
    if x < 0: return -1
    return 0

class Plane:
    """A plane in cartesian 3D space"""
    def __init__(self, a, b, c, d):
        self.a = float(a)
        self.b = float(b)
        self.c = float(c)
        self.d = float(d)

    def dot(self, p):
        x, y, z = map(float, p)
        return self.a*x + self.b*y + self.c*z + self.d

    def intersection(self, p1, p2):
        """Returns the coordinates of intersection of the line formed by points
        p1 and p2 with the plane."""
        x1, y1, z1 = map(float, p1)
        x2, y2, z2 = map(float, p2)
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1
        a, b, c, d = self.a, self.b, self.c, self.d
        den = (a*dx + b*dy + c*dz)
        if den == 0:
            raise RuntimeError("Could not compute intersection. Line appears to be parallel to the plane.")
        t = - (a*x1 + b*y1 + c*z1 + d)/den
        if t < 0 or t > 1:
            raise RuntimeError("Intersection does not fall within.")
        x = x1 + t * dx
        y = y1 + t * dy
        z = z1 + t * dz
        #print("Intersection:", dx, dy, dz)
        return (x, y, z)

class PSection(Plane):
    def compute_crossings(self, trajectory):
        top2bot = []
        bot2top = []
        x, y, z = trajectory.x, trajectory.y, trajectory.z
        N = len(x)
        p_prev = (x[0], y[0], z[0])
        s_prev = sign(self.dot(p_prev))
        for i in range(1, N):
            p = x[i], y[i], z[i]
            s = sign(self.dot(p))
            #print(s_prev, s, p)
            if s_prev != s:
                #print("Intersection within:", p_prev, p)
                p_int = self.intersection(p_prev, p)
                if s_prev > 0:
                     top2bot.append(p_int)
                else:
                     bot2top.append(p_int)
            p_prev = p
            s_prev = s
        return top2bot, bot2top

def test_crossings():
    import math
    def traj_circle(radius, th_max = 10*math.pi, intervals = 100):
        import numpy as np
        import dynsystem
        th = np.linspace(0, th_max, intervals)
        s = radius * th
        x = radius * np.cos(th)
        y = radius * np.sin(th)
        z = np.zeros_like(th)
        return dynsystem.Trajectory(np.array([th, s, x, y, z]).transpose())
    p = PSection(2, 1, 4, 70)
    t = traj_circle(5.0, intervals = 5000)
    top2bot, bot2top = p.compute_crossings(t)
    print "Top to bottom:"
    for p in top2bot:
        print p
    print "Bottom to top:"
    for p in bot2top:
        print p
    import matplotlib.pyplot as plt
    plt.plot(t.x, t.y, 'ro-')
    plt.show()

if __name__ == '__main__':
    test_crossings()
