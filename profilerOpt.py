#!/usr/bin/env python3
#
# -*- coding: utf-8 -*-
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

import argparse
from sys import argv, stdout, stderr
from modules.readopts import preprocess_args, argparse2opts, readinput
from modules.opts import *
from modules.vbga import VBGA, DefaultSelectionMethod, RoulettSelectionMethod, RankSelectionMethod, TournamentSelectionMethod
from modules.stpParser import parseStpFile
from modules.multiprofile import multiProfile
from modules.configuration import ensemble
from modules.writetraj import *
from modules.weightcalculator import *
from datetime import datetime
from platform import platform, python_version
from os.path import isfile
from random import seed
import numpy as np

def printheader (fp):
    version = '1.0'
    prog    = 'profilerOpt'
    bars    = "-"*(len(prog)+len(version)+6)
    string = "                         *** {0} (v. {1}) ***\n".format(prog,version) 
    string+= "                             {0}\n\n".format(bars)
    string+= "Running in {0} with Python {1}.\n".format(platform(), python_version())
    string+= "Started execution at {0}.\n".format(datetime.strftime(datetime.now(),"%Y/%m/%d %H:%M:%S")) 
    string+= "Using {0} processors for parallelization.\n".format(cmdlineOpts.nProcs) 
    print(string, file=fp)

def printfooter (fp):
    string = "Finished execution at {0}.\n".format(datetime.strftime(datetime.now(),"%Y/%m/%d %H:%M:%S")) 
    print(string, file=fp)

def checkfiles (fnlist):
    for fn in fnlist:
        if not (isfile(fn)):
            raise IOError ("File {0} does not exist!".format(fn))

def selectionWrapper (seltype):
    if (seltype) == 1:
        return RoulettSelectionMethod
    elif (seltype) == 2:
        return RankSelectionMethod
    elif (seltype) == 3:
        return TournamentSelectionMethod
    else:
        raise RuntimeError ("Selection type {} not allowed.".format(seltype))

def crossWrapper (crosstype):
    # if   (crosstype) == 1:
    #     return multiProfile.uniformCrossover
    # elif (crosstype) == 2:
    #     return multiProfile.averageCrossover
    if (crosstype) == 1:
        return multiProfile.arithmeticCrossover
    elif (crosstype) == 2:
        return multiProfile.heuristicCrossover
    else:
        raise RuntimeError ("Cross-over type {} not allowed.".format(crosstype))

if __name__ == '__main__':

    progdescr = """
    profilerOpt is a Python program for simultaneous optimization of
    torsional and 1-4 Lennard-Jones parameters. It is based on simple
    in-house Python libraries for parameter optimization (using a
    genetic algorithm) and calculation of molecular mechanics energies
    and forces.
    """

    parser = argparse.ArgumentParser(description=progdescr, formatter_class=argparse.RawTextHelpFormatter)

    # hack to avoid printing "optional arguments:" in help message
    parser._optionals.title = "options"

    parser.add_argument('-np', dest='nProcs', required=False, type=int, default=1, help=
            "Number of subprocesses spawned in parallelization.")

    parser.add_argument('-r', metavar=('REF_1','REF_2'), dest='ref', nargs='+', required=True, type=str, help=
            "Reference data files.")

    parser.add_argument('-c', metavar=('COORDS_1','COORDS_2'),
                        dest='coords', nargs='+', required=True, type=str, help=
                "Torsional-scan trajectory files (.g96/.xyz/.gro).")

    parser.add_argument('-t', metavar=('STP_1', 'STP_2'), dest='pars', nargs='+', required=True, type=str, help=
                "Special-topology files (.stp).")

    parser.add_argument('-i', metavar='INP', dest='input', required=True, type=str, help=
                "Input file containing profilerOpt parameters (.inp).")

    parser.add_argument('-w', metavar=('WEI_1','WEI_2'), dest='wei', nargs='+', required=False, type=str, help=
            "Weight files (default = derive weights from INP file).")

    parser.add_argument('-op', metavar='PREFIX', dest='out', required=True, type=str, help=
                "Prefix for output files.")

    args = parser.parse_args(preprocess_args(argv[1:]))

    # map arguments to xxxOpts
    argparse2opts (args)

    # starting program
    printheader(stdout)

    # read parameters from input and put them in the xxxOpts classes
    readinput (cmdlineOpts.inputFile)

    # set seed
    if randOpts.seed == -1:
        seed()
    else:
        seed(randOpts.seed)

    # check if files exist or raise IOError
    checkfiles (cmdlineOpts.stpFiles + cmdlineOpts.refFiles + cmdlineOpts.trajFiles + cmdlineOpts.weiFiles)
    
    # read ref and wei data
    optOpts.refData = [np.loadtxt(r) for r in cmdlineOpts.refFiles]
    if (len(cmdlineOpts.weiFiles) > 0):
        optOpts.weiData = [np.loadtxt(w) for w in cmdlineOpts.weiFiles]
    else:
        # fill wei data with appropriate weights
        weiCalcObj = initializeWeightCalculator(optOpts.wTemp)
        weiCalcObj.setEnergiesFromFiles(cmdlineOpts.refFiles)
        optOpts.weiData = weiCalcObj.computeWeights()

    # read stp data into xxxOpts class
    optOpts.stpData = [parseStpFile(s, prepareOpt=True) for s in cmdlineOpts.stpFiles]

    # initialize writeTraj
    pars_fn = "{}.trp".format(cmdlineOpts.outPrefix)
    ene_fn  = "{}".format(cmdlineOpts.outPrefix)
    fit_fn = "{}.trr".format(cmdlineOpts.outPrefix)
    trajWriter = writer(pars_fn, ene_fn, fit_fn, writetrajOpts.nstp, writetrajOpts.nste)

    # initialize GA
    GA = VBGA(multiProfile, multiProfile.rmsdToData,
          selectionWrapper(vbgaOpts.selectType), vbgaOpts.selectNum, vbgaOpts.ntel,
          crossWrapper(vbgaOpts.crossType), vbgaOpts.crossRate,
          multiProfile.mutate, vbgaOpts.mutRate,
          vbgaOpts.popSize,
          cmdlineOpts.nProcs,
          trajWriter)

    # run GA
    GA.run(vbgaOpts.nGens)

    # write best ind
    optind = GA.getBest()
    optind.saveTraj(cmdlineOpts.outPrefix, 'xyz')
    optind.saveProfile(cmdlineOpts.outPrefix)
    optind.saveParameters(cmdlineOpts.outPrefix + '.ifp')
    
    # now, update minim, do final optimization and save
    print("\nPerforming final optimization... ", file=stdout, end='')
    optind.resetMinim(lastMinimOpts.minimType, lastMinimOpts.dx0, lastMinimOpts.dxm, lastMinimOpts.dele, lastMinimOpts.maxSteps)
    optind.minimizeProfiles(useWei=False)
    optind.saveTraj(cmdlineOpts.outPrefix + '_minim', 'xyz')
    optind.saveProfile(cmdlineOpts.outPrefix + '_minim')
    print("Done.", file=stdout)

    printfooter(stdout)
