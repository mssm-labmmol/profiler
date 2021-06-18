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

from .opts import optOpts, randOpts, minimOpts, lastMinimOpts, vbgaOpts, dihrestrOpts, geomcheckOpts, cmdlineOpts
from .randomizer import RandomizerFactory, LimiterDecorator, SignReverserDecorator
from sys import stderr
from .configuration import ensemble
import numpy as np
from functools import partial
import itertools

def getline (stream, comm = ';'):
    buff = ''
    while buff == '':
        buff = stream.readline()
        if buff == '':
            return buff
        if buff.strip() == '':
            return "END"
        buff = buff.strip()
        buff = buff.split(comm)[0]
    return buff

def getblock (stream, block_cont, comm = ';'):
    block_cont.clear()
    _buff = ''
    first = True
    while (True):
        _buff = getline(stream, comm)
        if (_buff != ''):
            if (_buff == 'END'):
                break
            block_cont += [fld for fld in _buff.split() if (fld != '[') and (fld != ']')]
            # check if all arguments are numbers
            offset = 3 if first else 0
            for x in _buff.split()[offset:]:
                try:
                    float(x)
                except ValueError:
                    raise ValueError("%s is not a numeric option for block %s." % (x, block_cont))
        else:
            return False
        first = False
    return True

def read_PARAMETEROPTIMIZATION (blockdict):
    blockname = 'parameter_optimization'
    blocklist = blockdict[blockname]
    optOpts.nTors = int(blocklist.pop(0))
    optOpts.nLJ   = int(blocklist.pop(0))
    optOpts.dihType = {1: 'standard', 2: 'ryckaert'}[int(blocklist.pop(0))]
    randOpts.mslots = True
    optOpts.optTors = []
    for j in range(optOpts.nTors):
        currKMask = [0, 0, 0, 0, 0, 0]
        currPhiMask = [0, 0, 0, 0, 0, 0]
        optTorsMask = [0, 0, 0, 0, 0, 0]
        for i in range(6):
            nt = int(blocklist.pop(0))
            if (nt == 0):
                pass
            elif (nt == 1):
                currKMask[i] = 1
                optTorsMask[i] = 1
            elif (nt == 2):
                currKMask[i] = 1
                currPhiMask[i] = 1
                optTorsMask[i] = 1
            else:
                raise Exception()
        optOpts.kMask.append(currKMask)
        optOpts.phiMask.append(currPhiMask)
        optOpts.optTors.append(optTorsMask)
        
    if optOpts.nTors != len(optOpts.kMask):
        raise ValueError("Input file: Number of dihedral types do not match specification.")

    optOpts.optTors = np.array(optOpts.optTors, dtype=np.uint8)

    for j in range(optOpts.nLJ):
        cs6switch = int(blocklist.pop(0))
        cs12switch = int(blocklist.pop(0))
        optOpts.LJMask.append((cs6switch, cs12switch))

    if optOpts.nLJ != len(optOpts.LJMask):
        raise ValueError("Input file: Number of LJ types do not match specification.")

    if (len(optOpts.LJMask) == 0):
        optOpts.LJMask = None
    if (len(optOpts.kMask) == 0):
        optOpts.kMask = None
    if (len(optOpts.phiMask) == 0):
        optOpts.phiMask = None
        
    optOpts.wTemp = float(blocklist.pop(0))

def read_PARAMETERRANDOMIZATION (blockdict):
    blockname = 'parameter_randomization'
    blocklist = blockdict[blockname]
    randOpts.torsPinv = int(blocklist.pop(0))
    if randOpts.torsPinv == 0:
        raise Exception("The value of PINV must be non-zero to allow not taking the phase constants into account; read the documentation for more details.")
    randOpts.torsDist = int(blocklist.pop(0))
    randOpts.torsMin  = float(blocklist.pop(0))
    randOpts.torsMax  = float(blocklist.pop(0))
    randOpts.torsMean  = float(blocklist.pop(0))
    randOpts.torsStddev  = float(blocklist.pop(0))
    randOpts.cs6Dist  = int(blocklist.pop(0))
    randOpts.cs6Min  = float(blocklist.pop(0))
    randOpts.cs6Max  = float(blocklist.pop(0))
    randOpts.cs6Mean  = float(blocklist.pop(0))
    randOpts.cs6Stddev  = float(blocklist.pop(0))
    randOpts.cs12Dist  = int(blocklist.pop(0))
    randOpts.cs12Min  = float(blocklist.pop(0))
    randOpts.cs12Max  = float(blocklist.pop(0))
    randOpts.cs12Mean  = float(blocklist.pop(0))
    randOpts.cs12Stddev  = float(blocklist.pop(0))

def read_SEED(blockdict):
    blockname = 'seed'
    blocklist = blockdict[blockname]
    randOpts.seed = int(blocklist.pop(0))

def read_EVOLUTIONARYSTRAT (blockdict):
    blockname = 'evolutionary_strat'
    blocklist = blockdict[blockname]
    strategy_switch = int(blocklist.pop(0))
    if strategy_switch == 1:
        vbgaOpts.strategy = 'CMA-ES'
    elif strategy_switch == 2:
        vbgaOpts.strategy = 'GA'
    else:
        raise ValueError
    vbgaOpts.popSize = int(blocklist.pop(0))
    if (vbgaOpts.popSize <= 0):
        raise Exception("Population size must be a positive number.")
    vbgaOpts.nGens = int(blocklist.pop(0))

def read_MINIMIZATION (blockdict):
    blockname = 'minimization'
    blocklist = blockdict[blockname]
    minimOpts.minimType = int(blocklist.pop(0))
    minimOpts.dx0 = float(blocklist.pop(0))
    minimOpts.dxm = float(blocklist.pop(0))
    minimOpts.maxSteps = int(blocklist.pop(0)) 
    minimOpts.dele = float(blocklist.pop(0))
    if (minimOpts.minimType != 0) and (minimOpts.maxSteps == 0):
        raise RuntimeError ("Minimization algorithm is non-zero but the number of steps is 0!")
    if (minimOpts.minimType == 0):
        minimOpts.minimType = 1
        minimOpts.maxSteps = 0
    lastMinimOpts.minimType = int(blocklist.pop(0))
    lastMinimOpts.dx0 = float(blocklist.pop(0))
    lastMinimOpts.dxm = float(blocklist.pop(0))
    lastMinimOpts.maxSteps = int(blocklist.pop(0)) 
    lastMinimOpts.dele = float(blocklist.pop(0))
    if (lastMinimOpts.minimType != 0) and (lastMinimOpts.maxSteps == 0):
        raise RuntimeError ("Minimization algorithm is non-zero but the number of steps is 0!")
    if (lastMinimOpts.minimType == 0):
        lastMinimOpts.minimType = 1
        lastMinimOpts.maxSteps = 0

def read_TORSIONALSCAN (blockdict):
    blockname = 'torsional_scan'
    blocklist = blockdict[blockname]
    dihrestrOpts.origin = int(blocklist.pop(0))
    dihrestrOpts.ntypes = int(blocklist.pop(0))
    if dihrestrOpts.ntypes < 1:
        raise RuntimeError("You need at least one dihedral-restraint type.")
    for type_ in range(dihrestrOpts.ntypes):
        dihrestrOpts.start.append(float(blocklist.pop(0)))
        dihrestrOpts.step.append( float(blocklist.pop(0)))
        dihrestrOpts.last.append(float(blocklist.pop(0)))
        dihrestrOpts.k.append(float(blocklist.pop(0)))
        dihrestrOpts.nPoints.append(int(int(dihrestrOpts.last[-1] - dihrestrOpts.start[-1])/dihrestrOpts.step[-1]) + 1)
        
def read_GEOMETRYCHECK (blockdict):
    blockname = 'geometry_check'
    blocklist = blockdict[blockname]
    geomcheckOpts.lvlBond = int(blocklist.pop(0))
    geomcheckOpts.tolBond = float(blocklist.pop(0))
    geomcheckOpts.lvlAng = int(blocklist.pop(0))
    geomcheckOpts.tolAng = float(blocklist.pop(0))
    geomcheckOpts.lvlRestr = int(blocklist.pop(0))
    geomcheckOpts.tolRestr = float(blocklist.pop(0))
    geomcheckOpts.lvlImpr = int(blocklist.pop(0))
    geomcheckOpts.tolImpr = float(blocklist.pop(0))

# returns a dictionary with the contents of the blocks with their names as keys
def input2dict (fn):
    fp = open(fn, 'r')
    out_dict = {}
    block_cont = []
    while (True):
        read_code = getblock(fp, block_cont)
        try:
            out_dict[block_cont[0]] = block_cont[1:]
        except IndexError:
            pass
        if not (read_code):
            break
    fp.close()
    return out_dict

def readinput (fn):
    outdict = input2dict(fn)
    read_PARAMETEROPTIMIZATION (outdict)
    read_PARAMETERRANDOMIZATION (outdict)
    read_EVOLUTIONARYSTRAT (outdict)
    read_MINIMIZATION (outdict)
    read_TORSIONALSCAN (outdict)
    read_SEED(outdict)

def preprocess_args (argv, file_flag = '-f'): # argv should not contain the script name
    first = True
    found_file_flag = False
    if (argv == []):
        print ("Argument list is empty!", file=stderr)
        exit()
    for i,arg in enumerate(argv):
        if (arg == file_flag):
            found_file_flag = True
            if not (first):
                print ("If present,", file_flag, "<file> must be the only option.", file=stderr)
                exit()
        else:
            if (found_file_flag):
                if (i != len(argv) - 1):
                    print ("If present,", file_flag, "<file> must be the only option.", file=stderr)
                    exit()
                # in this case all arguments passed to arparser are in the argfile
                args = []
                readargfile(arg, args)
                return args
        first = False
    return argv

def readargfile (fn, block_cont, comm = '#'):
    block_cont.clear()
    _buff = ''
    stream = open(fn, 'r')
    while (True):
        _buff = getline(stream, comm)
        if (_buff != ''):
            block_cont += _buff.split()
        else:
            break
    stream.close()

def argparse2opts (args):
    cmdlineOpts.nProcs = args.nProcs
    cmdlineOpts.trajFiles = args.coords
    cmdlineOpts.stpFiles = args.pars
    cmdlineOpts.weiFiles = args.wei if args.wei is not None else []
    cmdlineOpts.refFiles = args.ref
    cmdlineOpts.outPrefix = args.out
    cmdlineOpts.inputFile = args.input
    cmdlineOpts.debugEmm = args.emm
    cmdlineOpts.dihspecFiles = args.dih_spec
    # consistency check
    nref = len(args.ref)
    ncoords = len(args.coords)
    npars = len(args.pars)
    nwei = len(args.wei) if args.wei is not None else ncoords
    if (nref == ncoords) and (ncoords == npars) and (npars == nwei):
        optOpts.nSystems = nref
    else:
        raise RuntimeError("Number of reference files, coordinate files, parameter files and "
                "weight files (if supplied) do not match.")

# ====================================================================
# Factories/other utils for initializing objects based on input
# parameters.
# ====================================================================

def MaskLooper(loop_func, LJMasks=None, kMasks=None, phiMasks=None):
    if not (len(kMasks) == len(phiMasks)):
        raise ValueError()
    type_count = 0
    out = []
    types = []
    if (LJMasks is not None):
        for C6Switch, C12Switch in LJMasks:
            if C6Switch != 0:
                out.append(loop_func('c6'))
                types.append(type_count)
            if C12Switch != 0:
                out.append(loop_func('c12'))
                types.append(type_count)
            type_count += 1
    type_count_save = type_count
    if (kMasks is not None):
        for sArr in kMasks:
            for m,s in enumerate(sArr):
                if s != 0:
                    out.append(loop_func('k', m+1))
                    types.append(type_count)
            type_count += 1
    type_count = type_count_save
    if (phiMasks is not None):
        for sArr in phiMasks:
            for m,s in enumerate(sArr):
                if s != 0:
                    out.append(loop_func('phi', m+1))
                    types.append(type_count)
            type_count += 1
    return out, types

class IndexListCreator:
    def __init__(self, optType, optAtoms, optPairs, optDihs):
        self.optType = optType
        self.optAtoms = optAtoms
        self.optPairs = optPairs
        self.optDihs = optDihs
        self.lenLJ = self._getLenLJ()

    def _getAtomOrPair(self, i):
        if (self.optType == 'atom'):
            return self.optAtoms[i]
        elif (self.optType == 'pair'):
            return self.optPairs[i]
        else:
            raise ValueError

    def _getLenLJ(self):
        if (self.optType == 'atom'):
            return len(self.optAtoms)
        elif (self.optType == 'pair'):
            return len(self.optPairs)
        else:
            return 0
        
    def _getDih(self, i):
        return self.optDihs[i]

    def get(self, idx):
        if idx >= self.lenLJ:
            return self._getDih(idx - self.lenLJ)
        else:
            return self._getAtomOrPair(idx)

def RandomizerSwitch2Type(sw):
    switch2type = {
        1: 'uniform',
        2: 'lognormal',
        3: 'gaussian',
        4: 'uniform_dim',
    }
    return switch2type[sw]

def RandomizerSwitch2Parameters(sw, min_, max_, mean, stdev):
    switch2pars = {
        1: {'min': min_, 'max': max_},
        2: {'mean': mean, 'stdev': stdev},
        3: {'mean': mean, 'stdev': stdev},
        4: {'min': min_, 'max': max_},
    }
    return switch2pars[sw]

def _InitRandomizerCore(typestr, m=None):
    if (typestr == 'c6'):
        return LimiterDecorator(RandomizerFactory(RandomizerSwitch2Type(randOpts.cs6Dist), **RandomizerSwitch2Parameters(randOpts.cs6Dist, randOpts.cs6Min, randOpts.cs6Max, randOpts.cs6Mean, randOpts.cs6Stddev)),
                         randOpts.cs6Min, randOpts.cs6Max)
    elif (typestr == 'c12'):
        return LimiterDecorator(RandomizerFactory(RandomizerSwitch2Type(randOpts.cs12Dist), **RandomizerSwitch2Parameters(randOpts.cs12Dist, randOpts.cs12Min, randOpts.cs12Max, randOpts.cs12Mean, randOpts.cs12Stddev)),
                         randOpts.cs12Min, randOpts.cs12Max)
    elif (typestr == 'k'):
        return SignReverserDecorator(LimiterDecorator(RandomizerFactory(RandomizerSwitch2Type(randOpts.torsDist), **RandomizerSwitch2Parameters(randOpts.torsDist, randOpts.torsMin, randOpts.torsMax, randOpts.torsMean, randOpts.torsStddev)), randOpts.torsMin, randOpts.torsMax), randOpts.torsPinv)
    elif (typestr == 'phi'):
        return RandomizerFactory('uniform', low=-180.00, up=+180.00)
        
def InitRandomizers():
    LJMasks  = optOpts.LJMask
    kMasks   = optOpts.kMask
    phiMasks = optOpts.phiMask
    
    rds = MaskLooper(_InitRandomizerCore, LJMasks, kMasks, phiMasks)[0]
    return rds


# ================ Dihedral-Restraint Reference Values ===============

def genRefDihedrals_Systematic(stpData, scanFirst, scanStep, scanLast):
    _refPhi = []
    for type_, reflist in enumerate(stpData['refdihedrals']):
        for d in reflist:
            scanValues = [scanFirst[type_] + i * scanStep[type_] for i in
                          range(int((scanLast[type_] - scanFirst[type_])/scanStep[type_]) + 1)]
            _refPhi.append(scanValues)
    
    out = np.array(list(itertools.product(*_refPhi))).T
    if out.ndim == 1:
        out = out.reshape(1, -1)
    return out
    

def genRefDihedrals_Trajectory(stpData, trajFile):
    ens = ensemble()
    ens.readFromTrajectory(trajFile)
    refPhi = []
    for type_, reflist in enumerate(stpData['refdihedrals']):
        for d in reflist:
            zeroIndexed = [ax - 1 for ax in stpData['propers'][0][d][:4]]
            phis = []
            for conf in ens:
                phis.append(conf.getDihedral(*zeroIndexed))
            refPhi.append(phis)
    refPhi = np.array(refPhi)
    if refPhi.ndim == 1:
        refPhi = refPhi.reshape(1, -1)
    return refPhi

def genRefDihedrals_ExternalFile(dihedralSpecFile):
    out = (np.loadtxt(dihedralSpecFile)).T
    if out.ndim == 1:
        out = out.reshape(1, -1)
    return out

def genRefDihedrals(system):
    if (dihrestrOpts.origin == 1):
        return genRefDihedrals_Systematic(optOpts.stpData[system],
                                          dihrestrOpts.start,
                                          dihrestrOpts.step,
                                          dihrestrOpts.last)
    elif (dihrestrOpts.origin == 2):
        return genRefDihedrals_Trajectory(optOpts.stpData[system],
                                          cmdlineOpts.trajFiles[system])
    elif (dihrestrOpts.origin == 3):
        return genRefDihedrals_ExternalFile(cmdlineOpts.dihspecFiles[system])
    else:
        raise ValueError
