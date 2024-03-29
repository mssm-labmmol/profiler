#
# This file is part of the profilerTools suite (see
# https://github.com/mssm-labmmol/profiler).
#
# Copyright (c) 2020 mssm-labmmol
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np
from numpy.linalg import norm
from math import sqrt
from .fastmath import fastCross


def calcDistance(x1, y1, z1, x2, y2, z2):
    diff = np.array([x2, y2, z2]) - np.array([x1, y1, z1])
    out = sqrt(np.dot(diff, diff))
    return out


def calcDistance2(x1, y1, z1, x2, y2, z2):
    diff = np.array([x2, y2, z2]) - np.array([x1, y1, z1])
    out = np.dot(diff, diff)
    return out


def calcNorm(vec):
    return calcDistance(vec[0], vec[1], vec[2], 0, 0, 0)


def calcAngle(x1, y1, z1, x2, y2, z2, x3, y3, z3):
    v1 = np.array([x1, y1, z1]) - np.array([x2, y2, z2])
    v2 = np.array([x3, y3, z3]) - np.array([x2, y2, z2])
    out = np.degrees(
        np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))
    return out

def calcAngleRadians(x1, y1, z1, x2, y2, z2, x3, y3, z3):
    v1 = np.array([x1, y1, z1]) - np.array([x2, y2, z2])
    v2 = np.array([x3, y3, z3]) - np.array([x2, y2, z2])
    return np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))


def calcCosine(x1, y1, z1, x2, y2, z2, x3, y3, z3):
    v1 = np.array([x1, y1, z1]) - np.array([x2, y2, z2])
    v2 = np.array([x3, y3, z3]) - np.array([x2, y2, z2])
    out = np.dot(v1, v2) / (sqrt(np.dot(v1, v1)) * sqrt(np.dot(v2, v2)))
    return out


# From: https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
def calcDihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = np.array([x1, y1, z1])
    p1 = np.array([x2, y2, z2])
    p2 = np.array([x3, y3, z3])
    p3 = np.array([x4, y4, z4])
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= sqrt(np.dot(b1, b1))
    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(fastCross(np.array([
        b1,
    ]), np.array([
        v,
    ])), w)
    out = np.degrees(np.arctan2(y, x))
    return out

def calcDihedralRadians(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = np.array([x1, y1, z1])
    p1 = np.array([x2, y2, z2])
    p2 = np.array([x3, y3, z3])
    p3 = np.array([x4, y4, z4])
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= sqrt(np.dot(b1, b1))
    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(fastCross(np.array([
        b1,
    ]), np.array([
        v,
    ])), w)
    return np.arctan2(y, x)

# Same as dihedral -- a different function only for clarity.
def calcImproper(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4):
    return calcDihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)


def wrapAngle(angleInDegrees):
    angleInFirstQuadrant = angleInDegrees % 360
    if (angleInFirstQuadrant <= 180.0):
        return np.radians(angleInFirstQuadrant)
    else:
        return np.radians(angleInFirstQuadrant - 360)


def wrapAngleDegree(angleInDegrees):
    angleInFirstQuadrant = angleInDegrees % 360
    if (angleInFirstQuadrant <= 180.0):
        return angleInFirstQuadrant
    else:
        return angleInFirstQuadrant - 360


def wrapAngles(angleInDegrees):
    return np.array([wrapAngle(x) for x in angleInDegrees])


def wrapAnglesDegrees(angleInDegrees):
    return np.array([wrapAngleDegree(x) for x in angleInDegrees])


def getBonds(x, y, z, bonds):
    out = []
    for b in bonds:
        out.append (calcDistance (x[b[0]-1], y[b[0]-1], z[b[0]-1],\
                      x[b[1]-1], y[b[1]-1], z[b[1]-1]) )
    return np.array(out)


def getAngles(x, y, z, angles):
    out = []
    for a in angles:
        out.append (calcAngle (x[a[0]-1], y[a[0]-1], z[a[0]-1],\
                      x[a[1]-1], y[a[1]-1], z[a[1]-1],\
                      x[a[2]-1], y[a[2]-1], z[a[2]-1]))
    return np.array(out)


def getDihedrals(x, y, z, dihedrals):
    out = []
    for d in dihedrals:
        out.append (calcDihedral (x[d[0]-1], y[d[0]-1], z[d[0]-1],\
                      x[d[1]-1], y[d[1]-1], z[d[1]-1],\
                      x[d[2]-1], y[d[2]-1], z[d[2]-1],\
                      x[d[3]-1], y[d[3]-1], z[d[3]-1]))
    return np.array(out)


def getImpropers(x, y, z, impropers):
    out = []
    for d in impropers:
        out.append(calcImproper (x[d[0]-1], y[d[0]-1], z[d[0]-1],\
                      x[d[1]-1], y[d[1]-1], z[d[1]-1],\
                      x[d[2]-1], y[d[2]-1], z[d[2]-1],\
                      x[d[3]-1], y[d[3]-1], z[d[3]-1]))
    return np.array(out)
