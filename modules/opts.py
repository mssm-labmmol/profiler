class cmdlineOpts (object):    

    nProcs = 1 # Number of processors. 
    trajFiles = [] # Trajectory (or single-configuration) files. 
    stpFiles = [] # Special topology files. 
    weiFiles = [] # Weight files. 
    refFiles = [] # Reference-data files. 
    inputFile = "" # Input parameters.
    outPrefix = "profopt" # Prefix for output files. 

class optOpts (object):

    nTors = 1 # Controls optimization of dihedrals 
    optTors = [1, 2, 3, 4, 5, 6] # List of terms to be optimized
    dihType = 'standard' # Dihedral type ('standard' or 'ryckaert')
    wTemp = 0.0 # Temperature for Boltzmann weight calculation. 0 = infinity.
    nLJ = 1 # Controls optimization of Lennard-Jones parameters 
    nSystems = 1 # Number of systems.
    weiData = [] # Weight data.
    refData = [] # Reference data.
    stpData = [] # Stp data.

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

    popSize = -1 # Population size. 
    nGens = 50 # Number of generations. 
    selectType = 0 # Selection operator. 
    selectNum  = 10 # Number of selected individuals. 
    ntel = 4 # Number of elitized individuals.
    crossType = 0 # Crossover operator. 
    crossRate = 20 # Crossover rate. 
    mutRate = 5 # Mutation rate. 

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

    start = 0 # Starting angle. 
    step = 10 # Step. 
    last = 360 # Last angles. 
    nPoints = 37 # How many points in the grid.
    #width = 0 # Width of the flat part of the potential. 
    k = 5000 # Force constant. 

class geomcheckOpts (object):

    lvlBond = 1 # Check bonds level. 
    tolBond = 0.1 # Bond tolerance. 
    lvlAng = 1 # Check angles level. 
    tolAng = 0.5 # Angle tolerance. 
    lvlRestr = 2 # Check restraints level. 
    tolRestr = 0.1 # Restraints tolerance. 
    lvlImpr = 1 # Check improper level. 
    tolImpr = 0.2 # Improper tolerance. 
