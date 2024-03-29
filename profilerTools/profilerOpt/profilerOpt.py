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
import numpy as np
import pkg_resources

from datetime import datetime
from platform import platform, python_version
from os.path import isfile
from random import seed
from sys import argv, stdout, stderr
from textwrap import fill

from ..readopts import preprocess_args, argparse2opts, readinput
from ..opts import *
from ..stpParser import parseStpFile
from ..multiprofile import multiProfile
from ..configuration import ensemble
from ..writetraj import *
from ..weightcalculator import *
from ..evolstrat import *
from ..linear_regression import *


def my_fill(text):
    return fill(text, width=55)


def useEvolStrat():
    if vbgaOpts.strategy is not None:
        return True
    else:
        return False

def evolStratFactory(*args, **kargs):
    evolutionaryStrategy = vbgaOpts.strategy
    if (evolutionaryStrategy == "GA"):
        return DEAP_GA(*args, **kargs)
    if (evolutionaryStrategy == "CMA-ES"):
        return DEAP_CMAES(*args, **kargs)


def printheader(fp):
    version = pkg_resources.require("profilerOpt")[0].version
    prog = 'profilerOpt'
    bars = "-" * (len(prog) + len(version) + 6)
    string = "                         *** {0} (v. {1}) ***\n".format(
        prog, version)
    string += "                             {0}\n\n".format(bars)
    string += "Running in {0} with Python {1}.\n".format(
        platform(), python_version())
    string += "Started execution at {0}.\n".format(
        datetime.strftime(datetime.now(), "%Y/%m/%d %H:%M:%S"))
    string += "Using {0} processors for parallelization.\n".format(
        cmdlineOpts.nProcs)
    print(string, file=fp)


def printfooter(fp):
    string = "Finished execution at {0}.\n".format(
        datetime.strftime(datetime.now(), "%Y/%m/%d %H:%M:%S"))
    print(string, file=fp)


def checkfiles(fnlist):
    for fn in fnlist:
        if not (isfile(fn)):
            raise IOError("File {0} does not exist!".format(fn))


def selectionWrapper(seltype):
    if (seltype) == 1:
        return RoulettSelectionMethod
    elif (seltype) == 2:
        return RankSelectionMethod
    elif (seltype) == 3:
        return TournamentSelectionMethod
    else:
        raise RuntimeError("Selection type {} not allowed.".format(seltype))


def crossWrapper(crosstype):
    # if   (crosstype) == 1:
    #     return multiProfile.uniformCrossover
    # elif (crosstype) == 2:
    #     return multiProfile.averageCrossover
    if (crosstype) == 1:
        return multiProfile.arithmeticCrossover
    elif (crosstype) == 2:
        return multiProfile.heuristicCrossover
    else:
        raise RuntimeError("Cross-over type {} not allowed.".format(crosstype))


def check_STP_Input_Consistency():
    for i, stp in enumerate(optOpts.stpData):
        if optOpts.nTors != len(stp['optdihedrals']):
            raise ValueError(
                "Number of opt. dihedral types in file {} should be {}.".
                format(cmdlineOpts.stpFiles[i], optOpts.nTors))
        if stp['opttype'] == 'pair':
            comp_LJ = len(stp['optpairs'])
        elif stp['opttype'] == 'atom':
            comp_LJ = len(stp['optatoms'])
        elif stp['opttype'] is None:
            comp_LJ = 0
        else:
            raise ValueError("Unexpected branch.")
        if optOpts.nLJ != comp_LJ:
            raise ValueError(
                "Number of opt. LJ types in file {} should be {}.".format(
                    cmdlineOpts.stpFiles[i], optOpts.nLJ))


class ProfilerOptRunner:
    def __init__(self, args):
        progdescr = """
        profilerOpt is a Python program for simultaneous optimization of
        torsional and 1-4 Lennard-Jones parameters taking into account
        several systems simultaneously.
        
        For each system, the NDIHS reference dihedrals of its torsional
        scan are defined in its input STP file. These reference dihedral
        angles are distributed into NGROUPS reference-dihedral groups,
        each specified by an individual [ refdihedrals ] block. The
        dihedral angles in the same group share the same value of
        restraint force constant, defined in the INP file (option -i). The
        number of reference dihedrals NDIHS may may not be the same for
        all systems, but the number of reference-dihedral groups NGROUPS
        must be the same.
    
        There are three ways to specify the values of the torsional-scan
        angles for each system: (i) explicitly, via NSCAN x NDIHS matrixes
        supplied for each system in the -s option (NSCAN is the length of
        the torsional scan for the system); (ii) keeping the
        reference-dihedral-angle values of the input trajectories; (iii)
        performing a (multidimensional) systematic scan for each system.
        Choosing any of these approaches involves a particular setting of
        option values in the INP file, which might have to be chosen in
        accordance with the contents of the STP file of each molecule. We
        recommend reading Chapter 4 of the documentation to set these
        options appropriately.
    
      """

        parser = argparse.ArgumentParser(
            description=progdescr,
            formatter_class=argparse.RawTextHelpFormatter)

        # hack to avoid printing "optional arguments:" in help message
        parser._optionals.title = "options"

        parser.add_argument(
            '-np',
            dest='nProcs',
            required=False,
            type=int,
            default=1,
            help=my_fill("Number of subprocesses spawned in parallelization."))

        parser.add_argument('-r',
                            metavar=('REF_1', 'REF_2'),
                            dest='ref',
                            nargs='+',
                            required=True,
                            type=str,
                            help=my_fill("Reference data files."))

        parser.add_argument(
            '-c',
            metavar=('COORDS_1', 'COORDS_2'),
            dest='coords',
            nargs='+',
            required=True,
            type=str,
            help=my_fill("Torsional-scan trajectory files (.g96/.xyz/.gro)."))

        parser.add_argument('-t',
                            metavar=('STP_1', 'STP_2'),
                            dest='pars',
                            nargs='+',
                            required=True,
                            type=str,
                            help=my_fill("Special topology files (.stp)."))

        parser.add_argument(
            '-i',
            metavar='INP',
            dest='input',
            required=True,
            type=str,
            help=my_fill(
                "Input file containing profilerOpt parameters (.inp)."))

        parser.add_argument(
            '-w',
            metavar=('WEI_1', 'WEI_2'),
            dest='wei',
            nargs='+',
            required=False,
            type=str,
            help=my_fill(
                "Weight files (default = derive weights from INP file)."))

        parser.add_argument('-op',
                            metavar='PREFIX',
                            dest='out',
                            required=True,
                            type=str,
                            help=my_fill("Prefix for output files."))

        parser.add_argument('--debug-mm',
                            dest='emm',
                            default=False,
                            action='store_true',
                            help=argparse.SUPPRESS)

        parser.add_argument(
            '-s',
            metavar=('SPEC_1', 'SPEC_2'),
            dest='dih_spec',
            nargs='+',
            required=False,
            type=str,
            help=my_fill("Torsional-scan dihedral-angle files."))

        args = parser.parse_args(preprocess_args(args))

        # map arguments to xxxOpts
        argparse2opts(args)

        # starting program
        printheader(stdout)

        # read parameters from input and put them in the xxxOpts classes
        readinput(cmdlineOpts.inputFile)

        # set seed
        if randOpts.seed == -1:
            seed()
        else:
            seed(randOpts.seed)
            np.random.seed(randOpts.seed)

        # check if files exist or raise IOError
        checkfiles(cmdlineOpts.stpFiles + cmdlineOpts.refFiles +
                   cmdlineOpts.trajFiles + cmdlineOpts.weiFiles)

        # read ref and wei data
        optOpts.refData = [np.loadtxt(r) for r in cmdlineOpts.refFiles]
        if (len(cmdlineOpts.weiFiles) > 0):
            optOpts.weiData = [np.loadtxt(w) for w in cmdlineOpts.weiFiles]
        else:
            # fill wei data with appropriate weights
            weiCalcObj = initializeWeightCalculator(optOpts.wTemp)
            weiCalcObj.setEnergiesFromFiles(cmdlineOpts.refFiles)
            optOpts.weiData = weiCalcObj.computeWeights()
        # put ref data at zero average
        optOpts.refData = [r - np.mean(r) for r in optOpts.refData]

        # read stp data into xxxOpts class
        optOpts.stpData = [
            parseStpFile(s, prepareOpt=True) for s in cmdlineOpts.stpFiles
        ]

        # if no minimization is requested during the GA execution,
        # calculate the MM energies of nonoptimized terms and store them
        # globally
        if (minimOpts.maxSteps == 0):
            # Create a dummy individual.
            dummyIndividual = multiProfile()
            # Get nonOpt energies
            optOpts.emmData = dummyIndividual.getNonoptEnergy()

        # check consistency between number of dihedral/LJ types in stp
        # files and in the input file
        check_STP_Input_Consistency()

        if useEvolStrat():
            # initialize strategy
            evolStrat = evolStratFactory(popSize=vbgaOpts.popSize)
    
            # run
            evolStrat.run(vbgaOpts.nGens, nprocs=cmdlineOpts.nProcs)
    
            # This is for debugging purposes - also write E_MM energies.
            if (cmdlineOpts.debugEmm):
                for i, ind in enumerate(evolStrat.population):
                    nonOpt = ind.getNonoptEnergy()
                    for j in range(nonOpt.shape[0]):
                        np.savetxt(
                            cmdlineOpts.outPrefix + '_' + str(i + 1) + '_' +
                            str(j + 1) + '_mm.dat', nonOpt[j, :])
    
            # write best ind
            optind = evolStrat.getBest()
            optind.saveTraj(cmdlineOpts.outPrefix, 'xyz')
        else:
            # LLS-SC
            # ------
            # Initialize a multiProfile object
            mp = multiProfile()
            # Initialize the LLS driver
            lls_sc = LLS_SC(mp)
            # Run
            lls_sc.run(llsOpts.max_cycles, llsOpts.max_dpar,
                       self.get_target_data(),
                       wei=self.get_weis(),
                       reg_center=llsOpts.reg_center,
                       lamb=llsOpts.lamb)

            optind = mp

        optind.saveProfile(cmdlineOpts.outPrefix)
        optind.saveTraj(cmdlineOpts.outPrefix, 'xyz')

        optind.saveParameters(cmdlineOpts.outPrefix + '.ifp')

        # now, update minim, do final optimization and save
        print("\nPerforming final optimization... ", file=stdout, end='')
        optind.prepareMinim(lastMinimOpts.minimType, lastMinimOpts.dx0,
                            lastMinimOpts.dxm, lastMinimOpts.dele,
                            lastMinimOpts.maxSteps)
        optind.minimizeProfiles(useWei=False)
        optind.saveTraj(cmdlineOpts.outPrefix + '_minim', 'xyz')
        optind.saveProfile(cmdlineOpts.outPrefix + '_minim')
        print("Done.", file=stdout)

        # store best
        self.optind = optind

        # save full data if requested
        if args.emm:
            for ip, profile in enumerate(self.optind):
                profile.mmCalc.calcForEnsembleAndSaveToFile(
                    profile.ensemble,
                    "{}_full_{}.dat".format(args.out, ip+1),
                    saveOnlyTotal=False)

        printfooter(stdout)


    def get_optimal_parameters(self):
        return self.optind.getOptimizableParameters()


    def get_optimal_individual(self):
        return self.optind

    def get_target_data(self):
        return np.array(optOpts.refData).flatten()

    def get_weis(self):
        return np.array(optOpts.weiData).flatten()


def main():
    job = ProfilerOptRunner(argv[1:])


if __name__ == '__main__':
    main()
