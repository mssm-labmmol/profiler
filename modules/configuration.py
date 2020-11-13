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
import geometrypy as gp
from .G96parser import parseG96traj
from .XYZparser import parseXYZtraj
from .GROparser import parseGROtraj
from .coordParser import calcDistance, calcAngle, calcCosine, calcDihedral, calcImproper
from math import sqrt

class configuration (object):

    def __init__ (self, pos=None):
        if (pos is None):
            self.pos = np.array([])
        else:
            self.pos = np.array(pos)

    def copy (self):
        newPositions = [x.copy() for x in self.pos]
        return configuration(newPositions)

    def __getitem__ (self, i):
        return self.pos[i]

    def __setitem__ (self, i, value):
        self.pos[i] = value

    def __iadd__ (self, other):
        if type(other) is np.ndarray:
            self.pos += other
        return self

    def __repr__ (self):
        out = ""
        for at in self.pos:
            out += ("[ %.4f %.4f %.4f ]\n" % (at[0], at[1], at[2]))
        return out

    def __add__ (self, other):
        newconf = self.copy() 
        if type(other) is np.ndarray:
            for i in range(len(self.pos)):
                newconf[i] += other[i]

        return newconf

    def put (self, positions):
        for i in range(len(self.pos)):
            self.pos[i] = positions[i]

    def getDisplacement (self, i, j):
        return (self.pos[i,:] - self.pos[j,:])

    def getDisplacements (self, ilist, jlist):
        return gp.calculateDisplacements(self.pos, ilist.astype(np.int32), jlist.astype(np.int32))

    def getDistance (self, i, j):
        return calcDistance (self.pos[i][0], self.pos[i][1], self.pos[i][2],\
                self.pos[j][0], self.pos[j][1], self.pos[j][2])

    def getDistances (self, ilist, jlist):
        return gp.calculateDistances(self.pos, ilist.astype(np.int32), jlist.astype(np.int32))

    def getDistances2 (self, ilist, jlist):
        return gp.calculateDistances2(self.pos, ilist.astype(np.int32), jlist.astype(np.int32))

    def getAngle (self, i, j, k):
        return calcAngle (self.pos[i][0], self.pos[i][1], self.pos[i][2],\
                self.pos[j][0], self.pos[j][1], self.pos[j][2],\
                self.pos[k][0], self.pos[k][1], self.pos[k][2])

    def getAngles (self, ilist, jlist, klist):
        return np.degrees(gp.calculateAngles(self.pos, ilist.astype(np.int32), jlist.astype(np.int32), klist.astype(np.int32)))

    def getCosine (self, i, j, k):
        return calcCosine (self.pos[i][0], self.pos[i][1], self.pos[i][2],\
                self.pos[j][0], self.pos[j][1], self.pos[j][2],\
                self.pos[k][0], self.pos[k][1], self.pos[k][2])

    def getSines (self, ilist, jlist, klist):
        return gp.calculateSines(self.pos, ilist.astype(np.int32), jlist.astype(np.int32), klist.astype(np.int32))

    def getCosines (self, ilist, jlist, klist):
        return gp.calculateCosines(self.pos, ilist.astype(np.int32), jlist.astype(np.int32), klist.astype(np.int32))

    def getDihedral (self, i, j, k, l):
        return calcDihedral (self.pos[i][0], self.pos[i][1], self.pos[i][2],\
                self.pos[j][0], self.pos[j][1], self.pos[j][2],\
                self.pos[k][0], self.pos[k][1], self.pos[k][2],\
                self.pos[l][0], self.pos[l][1], self.pos[l][2])

    def getDihedrals (self, ilist, jlist, klist, llist):
        return np.degrees(gp.calculateDihedrals(self.pos, ilist.astype(np.int32), jlist.astype(np.int32), klist.astype(np.int32), llist.astype(np.int32)))

    def getImproper (self, i, j, k, l):
        return calcImproper (self.pos[i][0], self.pos[i][1], self.pos[i][2],\
                self.pos[j][0], self.pos[j][1], self.pos[j][2],\
                self.pos[k][0], self.pos[k][1], self.pos[k][2],\
                self.pos[l][0], self.pos[l][1], self.pos[l][2])

    def getImpropers (self, ilist, jlist, klist, llist):
        return self.getDihedrals(ilist.astype(np.int32), jlist.astype(np.int32), klist.astype(np.int32), llist.astype(np.int32))

    def shiftParticle (self, i, shiftVector):
        self.pos[i] += shiftVector

    def size(self):
        return len(self.pos)

    def writeToStream (self, fp, elements, fmt='xyz'):
        nparticles = len(self.pos)
        if (fmt == 'xyz'):
            fp.write("%d\n" % nparticles)
            fp.write("dihedral scan frame\n")
            for pos,el in zip(self.pos, elements):
                fp.write("%2s    %18.7f%18.7f%18.7f\n" % (el,10*pos[0],10*pos[1],10*pos[2]))
        elif (fmt == 'g96'):
            fp.write("TIMESTEP\n%15d%15.4f\nEND\n" % (0, 0))
            fp.write("POSITIONRED\n")
            for pos in self.pos:
                fp.write("%15.9f%15.9f%15.9f\n" % (pos[0],pos[1],pos[2]))
            fp.write("END\n")
            fp.write("BOX\n%15.9f%15.9f%15.9f\nEND\n" % (0,0,0))
        else:
            raise ValueError("Format {} not supported for output.".format(fmt))

    def writeToFile (self, fn):
        fp = open(fn, 'w')
        self.writeToStream(fp)
        fp.close()

class ensemble (object):

    def __init__ (self, initConfs=None, elements=None):
        if (initConfs is None):
            self.confs = []
        else:
            self.confs = initConfs
        if (elements is None):
            self.elements = []
        else:
            self.elements = elements

    def copy (self):
        newConfs = [x.copy() for x in self.confs]
        return ensemble(newConfs)

    def __getitem__ (self, i):
        return self.confs[i]

    def __setitem__ (self, i, value):
        self.confs[i] = value

    def __repr__ (self):
        out = ""
        for i,x in enumerate(self.confs):
            out += "member %d \n" % (i+1)
            out += x.__repr__()
        return out

    def appendConf (self, conf):
        self.confs.append( conf.copy() )

    def readFromTrajectory (self, trajFile):
        if (trajFile.endswith('.g96')) or (trajFile.endswith('.trc')):
            traj = parseG96traj(trajFile, [], [], [], [], [])[-1]
            elem = ['X'] * traj.shape[2]
        if (trajFile.endswith('.gro')):
            (traj, elem) = parseGROtraj(trajFile)
        if (trajFile.endswith('.xyz')):
            (traj, elem) = parseXYZtraj(trajFile)

        if elem is not None:
            self.elements = elem
        for conf in traj:
            xs = conf[0]
            ys = conf[1]
            zs = conf[2]
            positions = []
            for x,y,z in zip (xs,ys,zs):
                positions.append(np.array([x,y,z]))
            self.confs.append(configuration(np.array(positions)))

    def size (self):
        return len(self.confs)

    def writeToFile (self, fn, fmt='xyz'):
        fp = open(fn, 'w')
        if (fmt == 'g96'):
            fp.write("TITLE\ndihedral-scan trajectory\nEND\n")
        for conf in self.confs:
            conf.writeToStream(fp, self.elements, fmt)
        fp.close()


