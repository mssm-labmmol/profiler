import numpy as np
import re
from  .coordParser import *

def isEndOfBlock (line):
    return (re.match(r"^END", line))

# Advances stream until a /^stringId/ matches.
def G96goto (stream, stringId):
    for line in stream:
        if re.match(r"^" + re.escape(stringId) + r"", line):
            return

def G96advanceTimestep (stream):
    G96goto (stream, "POSITION")

def advanceAndReadCoords (stream):
    x = []
    y = []
    z = []
    G96advanceTimestep (stream)
    for line in stream:
        line = line.split('#')[0]
        if (line == ''):
            continue
        if (isEndOfBlock(line)):
            break
        else:
            [tx, ty, tz] = [float(v) for v in line.split()]
            x.append(tx)
            y.append(ty)
            z.append(tz)
    return (x,y,z)

# Returns a list [bondDistances, angleValues, dihedralValues, improperValues, distanceValues, positions]
# where e.g. angleValues[i][j] is the value of the angles 'j' in timestep 'i'.
#
def parseG96traj (g96filename, bonds, angles, dihedrals, impropers, nbPairs):
    positions = []
    bondDistances = []
    angleValues = []
    dihedralValues = []
    improperValues = []
    distanceValues = []
    fp = open(g96filename, "r")
    while (True):
        (x,y,z) = advanceAndReadCoords (fp)
        # This means that no Timesteps were found.
        if (x == []):
            break
        # else get bonds angles dihedrals impropers from x,y,z data
        # these functions come from coordParser.py
        bondDistances.append( getBonds (x, y, z, bonds) )
        angleValues.append( getAngles (x, y, z, angles) )
        dihedralValues.append(getDihedrals (x, y, z, dihedrals))
        improperValues.append(getImpropers (x, y, z, impropers))
        # getBonds can be used for NB also
        distanceValues.append( getBonds(x,y,z,nbPairs))
        positions.append( np.array([x,y,z]) )
    fp.close()
    basicOut = [bondDistances, angleValues, dihedralValues, improperValues, distanceValues, positions]
    return [np.array(v) for v in basicOut]
