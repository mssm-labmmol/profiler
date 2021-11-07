# This is a wrapper to replace geometrypy functions in case the build fails
# or is not desired.

try:
    from geometrypy import *
except ImportError:
    import numpy as np
    from math import sqrt
    from .slow_math_warning import slow_math_warning
    from .coordParser import calcAngleRadians, calcDihedralRadians

    slow_math_warning()

    # auxiliary
    def getAngles(pos, ilist, jlist, klist):
        n = len(ilist)
        out = np.zeros(n)
        for i in range(n):
            out[i] = getAngle(pos, ilist[i], jlist[i], klist[i])
        return out

    # auxiliary
    def getAngle(pos, i, j, k):
        return calcAngleRadians(pos[i][0], pos[i][1], pos[i][2], pos[j][0],
                                pos[j][1], pos[j][2], pos[k][0], pos[k][1],
                                pos[k][2])

    # auxiliary
    def getDihedral(pos, i, j, k, l):
        return calcDihedralRadians(pos[i][0], pos[i][1], pos[i][2], pos[j][0],
                                   pos[j][1], pos[j][2], pos[k][0], pos[k][1],
                                   pos[k][2], pos[l][0], pos[l][1], pos[l][2])


    # use numpy einsum
    einsum = np.einsum

    def crossProduct(v1, v2):
        dims = len(v1.shape)
        eijk = np.zeros((3, 3, 3))
        eijk[0, 1, 2] = eijk[1, 2, 0] = eijk[2, 0, 1] = 1
        eijk[0, 2, 1] = eijk[2, 1, 0] = eijk[1, 0, 2] = -1
        if dims == 2:
            return np.einsum('ijk,aj,ak->ai', eijk, v1, v2)
        elif dims == 1:
            out = np.einsum('ijk,aj,ak->ai', eijk, [v1], [v2])
            return out[0]
        else:
            raise Exception

    def calculateDisplacements(pos, ilist, jlist):
        n = len(ilist)
        return np.array([pos[ilist[i],:] - pos[jlist[i],:] for i in
                         range(n)])

    def calculateDistances(pos, ilist, jlist):
        n = len(ilist)
        out = np.zeros(n)
        for i in range(n):
            dij = pos[ilist[i],:] - pos[jlist[i],:]
            out[i] = sqrt(np.dot(dij, dij))
        return out

    def calculateDistances2(pos, ilist, jlist):
        n = len(ilist)
        out = np.zeros(n)
        for i in range(n):
            dij = pos[ilist[i]] - pos[jlist[i]]
            out[i] = np.dot(dij, dij)
        return out

    def calculateAngles(pos, ilist, jlist, klist):
        n = len(ilist)
        out = np.zeros(n)
        for i in range(n):
            out[i] = getAngle(pos, ilist[i], jlist[i], klist[i])
        return out

    def calculateSines(pos, ilist, jlist, klist):
        return np.sin(np.radians(getAngles(pos, ilist, jlist, klist)))

    def calculateCosines(pos, ilist, jlist, klist):
        return np.cos(np.radians(getAngles(pos, ilist, jlist, klist)))

    def calculateDihedrals(pos, ilist, jlist, klist, llist):
        n = len(ilist)
        out = np.zeros(n)
        for i in range(n):
            out[i] = getDihedral(pos, ilist[i], jlist[i], klist[i], llist[i])
        return out
