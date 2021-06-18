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

from .multiprofile import multiProfile
from .opts import optOpts
import numpy as np


class writer(object):
    def __init__(self, pars_fn, ene_prefix, fit_fn, pars_freq, ene_freq):
        """
        Parameters:
          pars_fn (string) - parameters-trajectory filename
          ene_fn (string) - energy-trajectory filename
          fit_fn (string) - fitness-trajectory filename
          pars_freq (int) - frequency to write to pars_fn
          ene_freq (int) - frequency to write to ene_fn
        """
        self.pars_fn = pars_fn
        self.fit_fn = fit_fn
        self.ene_prefix = ene_prefix
        pars_stream = open(pars_fn, 'w')
        pars_stream.close()
        fit_stream = open(fit_fn, 'w')
        fit_stream.write("# GEN               AVG              BEST\n")
        fit_stream.close()
        for i in range(optOpts.nSystems):
            ene_fn = "{}_{}.tre".format(self.ene_prefix, i + 1)
            ene_stream = open(ene_fn, 'w')
            ene_stream.close()
        self.pars_freq = pars_freq
        self.ene_freq = ene_freq
        # just to create the file!

    def write_to_fit(self, generation, avg_fit, best_fit):
        """
        Parameters:
          generation (int) - generation we're in
          avg_fit (float) - average fitness of the population
          best_fit (float) best fitness of the population
        """
        fit_stream = open(self.fit_fn, 'a')
        fit_stream.write("{:>5d}{:>18.7e}{:>18.7e}\n".format(
            generation + 1, avg_fit, best_fit))
        fit_stream.close()

    def write_to_ene(self, generation, population):
        """
        Parameters:
          generation (int) - generation we're in
          population (list of multiProfile) - population of multiProfile individuals
        """
        for k in range(optOpts.nSystems):
            ene_fn = "{}_{}.tre".format(self.ene_prefix, k + 1)
            ene_stream = open(ene_fn, 'a')
            ene_stream.write("{:<5d}\n".format(generation + 1))
            ene_stream.write("# {:>6}".format("ANG"))
            for i in range(len(population)):
                ind_str = "IND{}".format(i + 1)
                ene_stream.write("{:>18}".format(ind_str))
            ene_stream.write("\n")
            refphi = population[0].profiles[k].refPhi
            for j in range(len(refphi)):
                ene_stream.write("{:>8.2f}".format(refphi[j]))
                for i in range(len(population)):
                    ene_stream.write("{:>18.7e}".format(
                        population[i].profiles[k].enerProfile[j]))
                ene_stream.write("\n")
            ene_stream.close()

    def write_to_pars(self, generation, population):
        """
        Parameters:
          generation (int) - generation we're in
          population (list of multiProfile) - population of multiProfile individuals
        """
        dihtype = optOpts.dihType
        pars_stream = open(self.pars_fn, 'a')
        pars_stream.write("{:<5d}\n".format(generation + 1))
        # torsions
        if (dihtype == 'standard'):
            for i in range(6):
                # pre-format
                k_str = "K{}".format(i + 1)
                phi_str = "PHI{}".format(i + 1)
                m_str = "M{}".format(i + 1)
                if (i == 0):
                    pars_stream.write("# {:>6}{:>18}{:>5}".format(
                        phi_str, k_str, m_str))
                else:
                    pars_stream.write("{:>8}{:>18}{:>5}".format(
                        phi_str, k_str, m_str))
        elif (dihtype == 'ryckaert'):
            for j in range(6):
                # pre-format
                c_str = "C{}".format(j)
                if (j == 0):
                    pars_stream.write("# {:>16}".format(c_str))
                else:
                    pars_stream.write("{:>18}".format(c_str))
        # lj
        if (optOpts.nLJ == -2):
            pars_stream.write("{:>18}".format("CS12"))
        if (optOpts.nLJ == -1):
            pars_stream.write("{:>18}".format("CS6"))
        if (optOpts.nLJ == 0):
            pass
        if (optOpts.nLJ == 1):
            pars_stream.write("{:>18}{:>18}".format("CS6", "CS12"))

        pars_stream.write("{:>18}\n".format("WRMSD"))
        for i in range(len(population)):
            if (dihtype == 'standard'):
                for j in range(6):
                    if (j in optOpts.optTors):
                        k = np.argwhere(optOpts.optTors == j)[0][0]
                        pars_stream.write("{:>8.2f}{:>18.7f}{:>5d}".format(
                            population[i].phi[k], population[i].k[k],
                            int(population[i].m[k])))
                    else:
                        pars_stream.write("{:>8.2f}{:>18.7f}{:>5d}".format(
                            0.0, 0.0, j + 1))
            elif (dihtype == 'ryckaert'):
                for j in range(6):
                    if (j in optOpts.optTors):
                        pars_stream.write("{:>18.7e}".format(
                            population[i].k[np.argwhere(
                                optOpts.optTors == j)[0][0]]))
                    else:
                        pars_stream.write("{:>18.7e}".format(0.0))
            if (optOpts.nLJ == -2):
                pars_stream.write("{:>18.7e}".format(population[i].cs12))
            if (optOpts.nLJ == -1):
                pars_stream.write("{:>18.7e}".format(population[i].cs6))
            if (optOpts.nLJ == 0):
                pass
            if (optOpts.nLJ == 1):
                pars_stream.write("{:>18.7e}{:>18.7e}".format(
                    population[i].cs6, population[i].cs12))
            pars_stream.write("{:>18.7e}\n".format(population[i].fitValue))
        pars_stream.close()

    def write(self, generation, population, avgFit, bestFit):
        self.write_to_fit(generation, avgFit, bestFit)
        if (self.pars_freq != 0) and (generation % self.pars_freq) == 0:
            self.write_to_pars(generation, population)
        if (self.ene_freq != 0) and (generation % self.ene_freq) == 0:
            self.write_to_ene(generation, population)
