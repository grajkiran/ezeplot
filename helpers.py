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

def is_inside(p, limits, strict = True):
    """If strict is True, then this function returns True only if all the
    limits are satisfied. Otherwise it returns True if any one limits is
    satisfied"""
    n_outside = 0
    if limits is None:
        return True
    for i in range(len(p)):
        if limits[i] is None:
            continue
        lower = min(limits[i])
        upper = max(limits[i])
        if p[i] < lower or p[i] > upper:
            n_outside += 1
    if strict:
        return n_outside == 0
    else:
        return n_outside < len(p)

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

def parse_coords(string):
    coords = list(map(float, string.replace(",", " ").split()))
    if len(coords) > 3:
        raise ValueError("Too many values.")
    coords.extend([0.0, 0.0, 0.0])
    return coords[:3]
