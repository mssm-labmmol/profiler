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


def parseGROtraj(fn):
    fp = open(fn, 'r')
    traj = []
    # file loop
    while (True):
        line = fp.readline()  # title - ignored
        if (line == ''):
            break
        line = fp.readline()  # number of atoms
        fields = line.split()
        natoms = int(fields[0])
        xs = []
        ys = []
        zs = []
        els = []
        for i in range(natoms):
            line = fp.readline()
            # gro is fixed-format
            ele = line[10:15].strip()[0]
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            els.append(ele)
            xs.append(x)
            ys.append(y)
            zs.append(z)
        line = fp.readline()  # box
        traj.append(np.array([xs, ys, zs]))
    fp.close()
    return (traj, els)
