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

class cmdlineOpts (object):    

    nProcs = 1 # Number of processors. 
    trajFiles = [] # Trajectory (or single-configuration) files. 
    stpFiles = [] # Special topology files. 
    weiFiles = [] # Weight files. 
    refFiles = [] # Reference-data files. 
    inputFile = "" # Input parameters.
    outPrefix = "profopt" # Prefix for output files. 
    dihspecFiles = None
    debugEmm = False

class optOpts (object):

    nTors = 1  # Controls optimization of dihedrals ; >0 is the number of torsional types.
    bOptPhase = False # Controls optimization of phases when dihType == 'standard'
    optTors = [1, 2, 3, 4, 5, 6] # List of terms to be optimized
    LJMask = [] # Mask of LJ terms to be optimized.
    kMask = [] # Mask of force constants to be optimized.
    phiMask = [] # Mask of phases to be optimized.
    dihType = 'standard' # Dihedral type ('standard' or 'ryckaert')
    wTemp = 0.0 # Temperature for Boltzmann weight calculation. 0 = infinity.
    nLJ = 1  # Controls optimization of Lennard-Jones parameters; >0 is the number of LJ types.
    nSystems = 1 # Number of systems.
    weiData = [] # Weight data.
    refData = [] # Reference data.
    stpData = [] # Stp data.
    emmData = None
    indexList = []  # 

class randOpts (object):

    seed = -1 # Seed
    jumpahead_n = 0 # Global jumpahead value.
    jumpahead_skip = 317 # Global jumpahead skip.
    torsDist = 2 # Distribution for randomization of dihedral force constants. 
    torsMin = 0 # Minimum value for dihedral force constants. 
    torsMax = 15 # Maximum value for dihedral force constants. 
    torsMean = 5.0 # Mean of distribution of dihedral force constant. 
    torsStddev = 15.0 # Standard deviation of distribution of dihedral force constant. 
    torsPinv = 50 # Probability of sign inversion of dihedral force constant. 
    cs6Dist = 2 # Distribution for randomization of CS6. 
    cs6Min = 0 # Minimimum value for CS6. 
    cs6Max = 1.0 # Maximimum value for CS6. 
    cs6Mean = 5.0 # Mean of distribution of CS6. 
    cs6Stddev = 15.0 # Standard deviation of distribution of CS6. 
    cs12Dist = 2 # Distribution for randomization of CS12. 
    cs12Min = 0 # Minimum value for CS12. 
    cs12Max = 1.0 # Maximum value for CS12. 
    cs12Mean = 5.0 # Mean of distribution of CS12. 
    cs12Stddev = 15.0 # Standard deviation of distribution of CS12. 
    mslots = True # Control of individual type.
    mslots_phi = 0.00 # Base phi value when mslots is used.

class vbgaOpts (object):

    strategy = "CMA-ES" # Method 
    popSize = -1 # Population size. 
    nGens = 50 # Number of generations. 

class writetrajOpts (object):

    nste = 5 # Freq to write profiles. 
    nstp = 5 # Freq to write parameters. 

class minimOpts (object):

    minimType = 1 # Minimization algorithm 
    dx0 = 0.01 # DX0 value. 
    dxm = 0.10 # Max DX. 
    maxSteps = 100 # Maximum number of iterations. 
    dele = 1.0e-09 # Energy threshold for convergence.

class lastMinimOpts (object):
    minimType = 1 # Minimization algorithm 
    dx0 = 0.01 # DX0 value. 
    dxm = 0.10 # Max DX. 
    maxSteps = 5000 # Maximum number of iterations. 
    dele = 1.0e-06 # Energy threshold for convergence.

class dihrestrOpts (object):

    origin = 1 # Origin of reference dihedral-angle values:
               # 1: Systematic
               # 2: From values in trajectory
               # 3: From values specified in external files
    ntypes = 0 # Number of dihedral-restraint types.
    start = [] # Starting angles.
    step = [] # Step angles.
    last = [] # Last angles. 
    nPoints = [] # How many points in each grid.
    #width = 0 # Width of the flat part of the potential. 
    k = [] # Force constants

class geomcheckOpts (object):

    lvlBond = 1 # Check bonds level. 
    tolBond = 0.1 # Bond tolerance. 
    lvlAng = 1 # Check angles level. 
    tolAng = 0.5 # Angle tolerance. 
    lvlRestr = 2 # Check restraints level. 
    tolRestr = 0.1 # Restraints tolerance. 
    lvlImpr = 1 # Check improper level. 
    tolImpr = 0.2 # Improper tolerance. 
