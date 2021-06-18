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

from .configuration import *
from .energy_force_vectorized import *
from .stpParser import *
from math import sqrt


class steepestDescentsMinimizer(object):
    def __init__(self, dx0, dxm, maxSteps, dele, mmCalc):
        self.dx0 = dx0
        self.dxm = dxm
        self.maxSteps = maxSteps
        self.dele = dele
        self.mmCalc = mmCalc

        self.energies = []
        self.ens = ensemble([])

    def clear(self):
        self.energies = []
        self.ens = ensemble([])
        self.mmCalc.clearForces()

    def appendEneConf(self, ene, conf):
        self.energies.append(ene)
        self.ens.appendConf(conf)

    def run(self, conf, debug=False):
        self.clear()
        dx = self.dx0
        dele = self.dele + 10
        iteration = 0
        (curr_ene,
         curr_f) = self.mmCalc.calcForConf(conf,
                                           calcForce=(self.maxSteps != 0),
                                           minim=True)
        self.appendEneConf(curr_ene, conf)
        while (True):
            if (iteration == self.maxSteps):
                break
            # dele is positive
            if (np.abs(dele) < self.dele):
                break
            (old_ene, old_f) = (curr_ene, curr_f)
            self.mmCalc.setForces(old_f)
            # Applicate the forces.
            self.mmCalc.normalizeForces()
            self.mmCalc.applyForcesToConfWithFactor(conf, dx)
            # Get new forces and energies.
            (curr_ene, curr_f) = self.mmCalc.calcForConf(conf, minim=True)
            # Calculate DELE.
            dele = curr_ene['total'] - old_ene['total']
            # If energy has increased, reject new configuration and
            # halve dx for the next step.
            if (dele > 0):
                # This effectively restores the old configuration.
                self.mmCalc.setForces(old_f)
                self.mmCalc.normalizeForces()
                self.mmCalc.applyForcesToConfWithFactor(conf, -dx)
                dx *= 0.50
                # This resets the new energies and forces to their old values.
                curr_ene = old_ene
                curr_f = old_f
                dele = self.dele + 10  # This is to avoid considering
                                       # the previous DELE in
                                       # convergence check.
            else:
                # Store new energies and configurations.
                self.appendEneConf(curr_ene, conf)
                # Increase DX, but keep it below DXM.
                dx *= 1.20
                if (dx > self.dxm):
                    dx = self.dxm
            iteration += 1

    def getUnrestrainedEnergy(self):
        return self.energies[-1]['total'] - self.energies[-1]['restraints']

    def getRestraintEnergy(self):
        return self.energies[-1]['restraints']

    def getFin(self):
        return self.ens[-1]

    def writeEnergiesToFile(self, fn, removeRestraints=False):
        if removeRestraints:
            np.savetxt(fn,
                       [x['total'] - x['restraints'] for x in self.energies])
        else:
            np.savetxt(fn, [x['total'] for x in self.energies])

    def writeEnsembleToTraj(self, fn):
        self.ens.writeToFile(fn)

    def runAndSave(self,
                   conf,
                   elements,
                   fn_ener,
                   fn_traj,
                   removeRestraints=False,
                   debug=False):
        self.run(conf, debug)
        self.writeEnergiesToFile(fn_ener, removeRestraints)
        self.ens.elements = elements
        self.writeEnsembleToTraj(fn_traj)


class conjugateGradientMinimizer(object):
    def __init__(self,
                 calc,
                 dx0=0.01,
                 dxm=0.10,
                 nsteps=100,
                 bfac=2.0,
                 prec=1.0e-08):
        self.dx0 = dx0
        self.dxm = dxm
        self.nsteps = nsteps
        self.bfac = bfac
        self.prec = prec
        self.energies = []
        self.ensemble = ensemble([])
        self.mmCalc = calc

    def clear(self):
        self.energies = []
        self.ensemble = ensemble([])
        self.mmCalc.clearForces()

    def run(self, conf):
        self.clear()
        stepEnergy, stepForce = self.mmCalc.calcForConf(conf, minim=True)
        # Store.
        self.ensemble.appendConf(conf)
        self.energies.append(stepEnergy)
        stepEnergy = stepEnergy['total']
        stepDir = stepForce

        # Outer cycle - actual number of steps.
        for step in range(self.nsteps):

            # Restore some quantities.
            dx = self.dx0
            ncyc = 0
            energyA = stepEnergy
            forceA = stepForce
            boundA = 0
            gA = np.sum(stepDir * forceA)
            boundB = dx / sqrt(np.sum(stepDir * stepDir))
            # Inner cycle - setting bounds a, b, including check on sMin.

            bcalcB = True
            while True:

                self.mmCalc.applyForcesToConfWithFactor(
                    conf, boundB - boundA, stepDir)

                if bcalcB:
                    energyB, forceB = self.mmCalc.calcForConf(conf, minim=True)
                    energyB = energyB['total']
                    gB = np.sum(stepDir * forceB)
                # Terminating conditions.
                if (energyB >= energyA) or (gB < 0):
                    # try interpolation
                    # Compute sMin.
                    try:
                        Z = 3 * (energyA - energyB) / (boundB -
                                                       boundA) - gA - gB
                        W = sqrt(Z * Z - gA * gB)
                        sMin = boundB - (boundB - boundA) * (W - Z - gB) / (
                            gA - gB + 2 * W)
                    except ValueError:
                        raise ValueError(
                            "Unexpected values were found during"
                            " conjugate-gradients minimization. Please, "
                            "try again with a smaller DX0."
                        )

                    # Energies and forces at sMin.
                    self.mmCalc.applyForcesToConfWithFactor(
                        conf, sMin - boundB, stepDir)
                    energySmin, forceSmin = self.mmCalc.calcForConf(conf,
                                                                    minim=True)
                    gMin = np.sum(forceSmin * stepDir)

                    if (energySmin['total'] <=
                            energyA) and (energySmin['total'] <= energyB):
                        break
                    else:
                        if gMin < 0:
                            # repeat for (a, sMin)
                            # move configuration to A
                            self.mmCalc.applyForcesToConfWithFactor(
                                conf, boundA - sMin, stepDir)
                            energyB = energySmin['total']
                            forceB = forceSmin.copy()
                            gB = gMin
                            boundB = sMin
                            ncyc += 1
                            # no need to calculate B
                            bcalcB = False
                        else:
                            # repeat for (sMin, b)
                            energyA = energySmin['total']
                            forceA = forceSmin.copy()
                            gA = gMin
                            boundA = sMin
                            ncyc += 1
                            # no need to calculate B
                            bcalcB = False
                else:
                    dx *= self.bfac
                    if (dx >= self.dxm):
                        raise Exception(
                            "DXM reached in conjugate-gradients method.")
                    # Save for next step.
                    energyA = energyB
                    forceA = forceB.copy()
                    gA = gB
                    boundA = boundB
                    boundB *= self.bfac
                    bcalcB = True
                    ncyc += 1

            # Store.
            self.ensemble.appendConf(conf)
            self.energies.append(energySmin)
            energySmin = energySmin['total']

            if (np.abs(energySmin - stepEnergy) <= self.prec):
                return

            # Save for next step.
            stepDir = forceSmin + stepDir * (np.sum(forceSmin * forceSmin) /
                                             np.sum(stepForce * stepForce))
            stepEnergy = energySmin
            stepForce = forceSmin.copy()
            prevEnergy = stepEnergy
            prevForce = stepForce.copy()

    def getUnrestrainedEnergy(self):
        return self.energies[-1]['total'] - self.energies[-1]['restraints']

    def getRestraintEnergy(self):
        return self.energies[-1]['restraints']

    def getFin(self):
        return self.ensemble[-1]

    def writeEnergiesToFile(self, fn, removeRestraints=False):
        if removeRestraints:
            np.savetxt(fn,
                       [x['total'] - x['restraints'] for x in self.energies])
        else:
            np.savetxt(fn, [x['total'] for x in self.energies])

    def writeEnsembleToTraj(self, fn):
        self.ensemble.writeToFile(fn)

    def runAndSave(self,
                   conf,
                   elements,
                   fn_ener,
                   fn_traj,
                   removeRestraints=False):
        self.run(conf)
        self.ensemble.elements = elements
        self.writeEnergiesToFile(fn_ener, removeRestraints)
        self.writeEnsembleToTraj(fn_traj)
