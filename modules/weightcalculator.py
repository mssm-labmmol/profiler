from abc import ABC, abstractmethod
import numpy as np

class weightStrategy(ABC):
    @abstractmethod
    def computeWeights(self, ener): pass

class uniformStrategy(weightStrategy):
    def __init__ (self): pass
    def computeWeights(self, ener):
        return [1.0 for e in ener]

class boltzmannStrategy(weightStrategy):
    def __init__(self, temp):
        self.temp = temp
        self.R    = 8.3145e-03
    def computeWeights(self, ener):
        min_ener = np.min(ener)
        return list(map(lambda x: np.exp(-(x - min_ener)/(self.R * self.temp)), ener))
        
class weightCalculator:

    def __init__(self, strategy):
        self.strategy = strategy

    def setEnergiesFromFiles(self, files):
        self.energies = [None for fn in files]
        for i,fn in enumerate(files):
            self.energies[i] = np.loadtxt(fn, usecols=(0,))

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
        raise Exception("Can't intialize weightCalculator from WGEN = {}".format(code))
    return weightCalculator(weiStrat)


        
