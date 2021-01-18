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

from .opts import *
from sys import stderr
import numpy as np

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
    for i in range(6):
        nt = int(blocklist.pop(0))
        if (nt == 0):
            pass
        elif (nt == 1):
            optOpts.optTors.append(i)
        else:
            raise Exception()
    optOpts.optTors = np.array(optOpts.optTors, dtype=np.uint8)
    optOpts.cs6switch = int(blocklist.pop(0))
    optOpts.cs12switch = int(blocklist.pop(0))
    if optOpts.nLJ == 1:
        if optOpts.cs6switch ==  1:
            if optOpts.cs12switch == 1:
                optOpts.nLJ = 1 # redundant, but yeah...
            else:
                optOpts.nLJ = -1
        else:
            if optOpts.cs12switch == 1:
                optOpts.nLJ = -2
            else:
                optOpts.nLJ = 0            
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

def read_WRITETRAJ (blockdict):
    blockname = 'writetraj'
    blocklist = blockdict[blockname]
    writetrajOpts.nste =   int(blocklist.pop(0))
    writetrajOpts.nstp =   int(blocklist.pop(0))

def read_GENETICALGORITHM (blockdict):
    blockname = 'genetic_algorithm'
    blocklist = blockdict[blockname]
    vbgaOpts.popSize = int(blocklist.pop(0))
    if (vbgaOpts.popSize <= 0):
        raise Exception("Population size must be a positive number.")
    vbgaOpts.nGens = int(blocklist.pop(0))
    vbgaOpts.selectType = int(blocklist.pop(0))
    vbgaOpts.ntel = int(blocklist.pop(0))
    vbgaOpts.crossType = int(blocklist.pop(0))
    vbgaOpts.crossRate = int(blocklist.pop(0))
    vbgaOpts.mutRate = int(blocklist.pop(0))
    vbgaOpts.selectNum  = vbgaOpts.popSize - vbgaOpts.ntel

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
    dihrestrOpts.start = float(blocklist.pop(0))
    dihrestrOpts.step  = float(blocklist.pop(0))
    dihrestrOpts.last  = float(blocklist.pop(0))
    dihrestrOpts.k     = float(blocklist.pop(0))
    dihrestrOpts.nPoints = int(int(dihrestrOpts.last - dihrestrOpts.start)/dihrestrOpts.step) + 1

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
        out_dict[block_cont[0]] = block_cont[1:]
        if not (read_code):
            break
    fp.close()
    return out_dict

def readinput (fn):
    outdict = input2dict(fn)
    read_PARAMETEROPTIMIZATION (outdict)
    read_PARAMETERRANDOMIZATION (outdict)
    read_WRITETRAJ (outdict)
    read_GENETICALGORITHM (outdict)
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
