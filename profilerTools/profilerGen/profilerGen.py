#!/usr/bin/env python3
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

from sys import argv, stdout, stderr
from datetime import datetime
from platform import platform, python_version
from os.path import isfile
from textwrap import fill

from ..readopts import (preprocess_args, genRefDihedrals_Systematic,
                        genRefDihedrals_Trajectory,
                        genRefDihedrals_ExternalFile)
from ..stpParser import parseStpFile
from ..configuration import *
from ..coordParser import calcDihedral, wrapAngleDegree
from ..energy_force_vectorized import MMCalculator
from ..minim import steepestDescentsMinimizer, conjugateGradientMinimizer
from ..multiprofile import profile as Profile


def my_fill(text):
    """Adjust text for a width of 55 characters (useful for argparse
    descriptions).

    """
    return fill(text, width=55)


def minimize_conf_with_profile(ensemble, out_prefix, stpData, refPhi,
                               restrConst, emAlgo, emDX0, emDXM, emDele,
                               emSteps):
    profile = Profile(stpData, None, refPhi, restrConst, emAlgo, emDX0, emDXM,
                      emDele, emSteps)
    if ensemble.size() <= 1:
        nconfs = profile.calculateNumberOfConfigurations()
        for i in range(nconfs - 1):
            ensemble.appendConf(ensemble[0])
    profile.replaceEnsemble(ensemble)
    profile.minimizeProfile(ensemble.elements)
    return profile


def printheader(fp):
    version = '1.0'
    prog = 'profilerGen'
    bars = "-" * (len(prog) + len(version) + 6)
    string = "                         *** {0} (v. {1}) ***\n".format(
        prog, version)
    string += "                             {0}\n\n".format(bars)
    string += "Running in {0} with Python {1}.\n".format(
        platform(), python_version())
    string += "Started execution at {0}.\n".format(
        datetime.strftime(datetime.now(), "%Y/%m/%d %H:%M:%S"))
    print(string, file=fp)


def printfooter(fp):
    string = "Finished execution at {0}.\n".format(
        datetime.strftime(datetime.now(), "%Y/%m/%d %H:%M:%S"))
    print(string, file=fp)


def checkfiles(fnlist):
    for fn in fnlist:
        if not (isfile(fn)):
            raise IOError("File {0} does not exist!".format(fn))


def debug_msg(msg, fp, verbose):
    if (verbose):
        print(msg, file=fp)


def err_msg(msg):
    raise RuntimeError(msg)


# Runner is a separate class just to facilitate testing Everything is
# done on initialization.  The only attribute is `self.energies',
# which is retrieved via `get_energies' for testing.
class ProfilerGenRunner:
    def __init__(self, args):
        progdescr = """
        profilerGen computes torsional-scan trajectories and their
        corresponding energy profiles.
    
        The NDIHS reference dihedrals of the torsional scan are defined in
        the input STP file. They are distributed into NGROUPS
        reference-dihedral groups, each specified by an individual
        [ refdihedrals ] block. The dihedral angles in the same group share
        the same value of restraint force constant, defined by the option
        -dk. This option must be supplied NGROUPS times, and it is applied
        to each reference-dihedral group following their order of
        appearance in the command string and in the STP file, respectively.
    
        There are three ways to specify the values of the torsional-scan
        angles. In the most flexible approach, the user supplies a NSCAN x
        NDIHS matrix file with the values of the reference dihedral angles
        for each torsional-scan configuration (NSCAN is the length of the
        torsional scan). This is done using the option -s. Alternatively,
        the user can supply parameters for a systematic multidimensional
        (dimension = NDIHS) scan with the option -dr. In this case, each
        -dr <RFRST> <RSTEP> <RLST> instance applies to all dihedral angles
        in the same reference-dihedral group, following their order of
        appearance in the command string and in the STP file,
        respectively. Each dihedral angle in a group is independently
        changed from RFRST to RLST (including RLST) in steps of RSTEP. The
        option -dr must be supplied NGROUPS times. If neither -s nor -dr
        is used, the torsional-scan angles are obtained from the input
        trajectory.
    
        """
        parser = argparse.ArgumentParser(
            description=progdescr,
            formatter_class=argparse.RawTextHelpFormatter)

        # hack to avoid printing "optional arguments:" in help message
        parser._optionals.title = "options"

        parser.add_argument(
            '-c',
            dest='coords',
            required=True,
            type=str,
            help=my_fill("Torsional-scan trajectory or single molecular "
                         "conformation (.xyz/.gro/.g96)."))

        parser.add_argument('-t',
                            metavar='STP',
                            dest='pars',
                            required=True,
                            type=str,
                            help=my_fill("Special topology file (.stp)."))

        parser.add_argument(
            '-dk',
            metavar=('RFCT_1', 'RFCT_2'),
            dest='dih_k',
            required=True,
            nargs='+',
            help=my_fill("Force constants for dihedral restraint (we suggest"
                         " 5000 kJ/(mol.rad^2)) for each type."))

        parser.add_argument('-op',
                            metavar='PREFIX',
                            dest='out',
                            required=True,
                            type=str,
                            help=my_fill("Prefix for output files."))

        parser.add_argument(
            '-s',
            metavar='SPEC',
            dest='dih_spec',
            type=str,
            required=False,
            help=my_fill("Torsional-scan angles specification file."))

        parser.add_argument(
            '-dr',
            metavar=('RFRST_1 RSTEP_1 RLST_1', 'RFRST_2 RSTEP_2 RLST_2'),
            action='append',
            dest='dih_restr',
            type=float,
            required=False,
            nargs='+',
            help=my_fill("Systematic-scan parameters for each type: "
                         "from RFRST to RLST with a RSTEP step."))

        parser.add_argument(
            '-min',
            metavar='MALG',
            dest='min_alg',
            required=False,
            type=int,
            default=1,
            help=my_fill(
                "Minimization algorithm: steepest descents (1 - default) "
                "or conjugate-gradient (2)."))

        parser.add_argument('-dx0',
                            metavar='DX0',
                            dest='min_dx0',
                            required=False,
                            type=float,
                            default=0.05,
                            help=my_fill(
                                "Initial step size DX0 in minimization"
                                " algorithm (default = 0.05 nm)."))

        parser.add_argument(
            '-dxm',
            metavar='DXM',
            dest='min_dxm',
            required=False,
            type=float,
            default=0.2,
            help=my_fill("Maximal step size DXM in"
                         " minimization algorithm (default = 0.20 nm)."))

        parser.add_argument('-dele',
                            metavar='DELE',
                            dest='min_dele',
                            required=False,
                            type=float,
                            default=1.0e-09,
                            help=my_fill(
                                "If |E_n+1 - E_n| <= DELE, minimization"
                                " stops (default = 1.0e-09)."))

        parser.add_argument(
            '-nsteps',
            metavar='NMAX',
            dest='min_nsteps',
            required=False,
            type=int,
            default=50000,
            help=my_fill(
                "Maximal number of steps in minimization (default = 50000)."))

        # this is in case arguments are provided in a -f <file>
        args = parser.parse_args(preprocess_args(args))

        # starting program
        printheader(stdout)
        # check if files exist or raise IOError
        checkfiles([args.coords, args.pars])

        if args.dih_restr is not None and (len(args.dih_restr) != len(
                args.dih_k)):
            raise RuntimeError(
                "Different number of types for -dr and -dk options.")

        restrConst = args.dih_k

        # get reference dihedrals from stp file
        stpdict = parseStpFile(args.pars)

        ens = ensemble([])
        ens.readFromTrajectory(args.coords)

        if args.dih_restr is None:
            if args.dih_spec is None:
                dihMode = 'trajectory'
            else:
                dihMode = 'external'
        else:
            if args.dih_spec is None:
                dihMode = 'systematic'
            else:
                raise Exception("Can't use -dr and -s at the same time.")

        if dihMode == 'systematic':
            scanFirst = [t_[0] for t_ in args.dih_restr]
            scanStep = [t_[1] for t_ in args.dih_restr]
            scanLast = [t_[2] for t_ in args.dih_restr]
            refPhi = genRefDihedrals_Systematic(stpdict, scanFirst, scanStep,
                                                scanLast)
            if (ens.size() > 1):
                debug_msg("Optimizing from input trajectory.", stdout, True)
            else:
                debug_msg("Optimizing from single configuration.", stdout,
                          True)
        elif dihMode == 'trajectory':
            refPhi = genRefDihedrals_Trajectory(stpdict, args.coords)
            debug_msg(
                "Optimizing from input trajectory (fixing reference dihedrals).",
                stdout, True)
        elif dihMode == 'external':
            refPhi = genRefDihedrals_ExternalFile(args.dih_spec)

        profile = minimize_conf_with_profile(ens, args.out, stpdict, refPhi,
                                             restrConst, args.min_alg,
                                             args.min_dx0, args.min_dxm,
                                             args.min_dele, args.min_nsteps)

        profile.ensemble.writeToFile(args.out + '.xyz', fmt='xyz')

        self.energies = profile.mmCalc.calcForEnsembleAndSaveToFile(
            profile.ensemble, "{}.dat".format(args.out))

        printfooter(stdout)

    def get_energies(self):
        return self.energies


def main():
    job = ProfilerGenRunner(argv[1:])


if __name__ == '__main__':
    main(argv[1:])
