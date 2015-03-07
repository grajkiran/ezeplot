#!/usr/bin/python3
import numpy as np
from math import sqrt

def scale_domain(d, factor, about = None):
    """Scales the domain by 'factor' maintaining 'about' at the same relative
    location from either ends. 'about' defaults to the midpoint of the interval.
    """
    lower = min(d)
    upper = max(d)
    if about is None:
        about = 0.5 * (lower + upper)
    delta_l = factor * (about - lower)
    delta_r = factor * (upper - about)
    return about - delta_l, about + delta_r

def is_inside(p, limits):
    if limits is None:
        return True
    for i in range(len(p)):
        if limits[i] is None:
            continue
        lower = min(limits[i])
        upper = max(limits[i])
        if p[i] < lower or p[i] > upper:
            return False
    return True

def freq_domain(times, amplitudes):
    """Perform and fft and return amplitude vs frequency."""
    from scipy.fftpack import fft
    dt = times[1]-times[0]
    T = dt*len(times)
    n_freqs = len(times)//2
    freqs = np.arange(n_freqs)/T
    amps = np.abs(fft(amplitudes)) / n_freqs
    return freqs, amps[:n_freqs]

def mag(vec):
    sum_sqr = 0.0
    for v in vec:
        sum_sqr += v*v
    return sqrt(sum_sqr)

def distance(p1, p2):
        sum_sqr = 0.0
        for i, j in zip(p1, p2):
            sum_sqr += (i-j)*(i-j)
        return sqrt(sum_sqr)

def curve_resample(x, y, count = 0, delta = 0.0, indices = False):
    from scipy.interpolate import interp1d
    if count <= 0 and delta <= 0.0:
        return np.arange(len(x), dtype = 'int'), x, y
    count_in = len(x)
    ind_in = np.arange(count_in)
    dist_in = np.zeros(count_in)
    for i in range(1, count_in):
        p1 = x[i-1], y[i-1]
        p2 = x[i], y[i]
        dist_in[i] = dist_in[i-1] + distance(p1, p2)
    if delta > 0.0:
        count_out = int(np.ceil(dist_in[-1]/delta))
    else:
        count_out = count
    if count_out <= 1:
        return np.arange(1), x[:1], y[:1]
    if count_out >= len(x):
        return np.arange(len(x), dtype = 'int'), x, y
    dist_out = np.linspace(0, dist_in[-1], count_out)
    ind_zero = interp1d(dist_in, ind_in, kind = 'zero')(dist_out)
    x_out = interp1d(dist_in, x)(dist_out)
    y_out = interp1d(dist_in, y)(dist_out)
    return np.array(ind_zero, dtype = 'int'), x_out, y_out

def curve_resample_old(points, count = 0, delta = 0.0):
    from scipy.interpolate import interp1d
    data = np.array(points)
    if count <= 0 and delta <= 0.0:
        return data
    count_in = len(data)
    ind_in = np.arange(count_in)
    x_in = data[:,0]
    y_in = data[:,1]
    dist_in = np.zeros(count_in)
    for i in range(1, count_in):
        dist_in[i] = dist_in[i-1] + distance(data[i], data[i-1])
    if delta > 0.0:
        count_out = int(np.ceil(dist_in[-1]/delta))
    else:
        count_out = count
    if count_out <= 1:
        return np.array(data[:1])
    dist_out = np.linspace(0, dist_in[-1], count_out)
    ind_out = interp1d(dist_in, ind_in)(dist_out)
    x_out = interp1d(ind_in, x_in)(ind_out)
    y_out = interp1d(ind_in, y_in)(ind_out)
    return np.transpose(np.array([x_out, y_out]))
