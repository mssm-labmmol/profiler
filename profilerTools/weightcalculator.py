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

from abc import ABC, abstractmethod
import numpy as np


class weightStrategy(ABC):
    @abstractmethod
    def computeWeights(self, ener):
        pass


class uniformStrategy(weightStrategy):
    def __init__(self):
        pass

    def computeWeights(self, ener):
        return [1.0 for e in ener]


class boltzmannStrategy(weightStrategy):
    def __init__(self, temp):
        self.temp = temp
        self.R = 8.3145e-03

    def computeWeights(self, ener):
        min_ener = np.min(ener)
        return list(
            map(lambda x: np.exp(-(x - min_ener) / (self.R * self.temp)),
                ener))


class weightCalculator:
    def __init__(self, strategy):
        self.strategy = strategy

    def setEnergiesFromFiles(self, files):
        self.energies = [None for fn in files]
        for i, fn in enumerate(files):
            self.energies[i] = np.loadtxt(fn, usecols=(0, ))

    def computeWeights(self):
        weis = []
        for ener in self.energies:
            weis.append(self.strategy.computeWeights(ener))
        return weis


def initializeWeightCalculator(code):
    if (code == 0.0):
        print("Using uniform weights.")
        weiStrat = uniformStrategy()
    elif (code > 0.0):
        print("Using Boltzmann weights with temp = {:f}".format(code))
        weiStrat = boltzmannStrategy(code)
    else:
        raise Exception(
            "Can't intialize weightCalculator from WGEN = {}".format(code))
    return weightCalculator(weiStrat)
