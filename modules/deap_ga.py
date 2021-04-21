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

from deap import base
from deap import creator
from deap import tools
from deap import algorithms

import sys

from .multiprofile import multiProfile
import multiprocessing
import numpy as np

def defaultInitializer():
    # Create first a multiprofile instance.
    mp = multiProfile()
    # Randomize it according to input.
    mp.randomize()
    # Return parameters.
    out = mp.getOptimizableParameters()
    return out

def defaultEvaluate(individual):
    mp = multiProfile()
    mp.setOptimizableParameters(slice(0, len(individual), 1), list(individual))
    mp.minimizeProfiles()
    out = mp.rmsdToData()
    return 1.0/out,

def defaultMut(individual, indpb=0.26, stdperc=0.1):
    """
    Each individual attribute suffers a Gaussian mutation with probability
    `indpb'. The standard deviations are `stdperc' times the corresponding mean.
    """ 
    mu = [0 for i in range(len(individual))]
    sigma = [stdperc * x for x in individual]
    return tools.mutGaussian(individual, mu, sigma, indpb)

def defaultCx(mother, father):
    return tools.cxOnePoint(mother, father)

def defaultSel(population, howMany):
    return tools.selRoulette(individuals=population,
                               k=howMany, 
                               fit_attr='fitness')

class DEAP_GA:

    def __init__(self,
                 popSize=200,
                 initFunc=defaultInitializer,
                 evalFunc=defaultEvaluate,
                 selFunc=defaultSel,
                 cxFunc=defaultCx,
                 mutFunc=defaultMut,
                 cxpb=0.23,
                 mutpb=0.39,
                 hofn=5,
                 nprocs=1):

        self.popSize = popSize
        self.cxpb = cxpb
        self.mutpb = mutpb

        toolbox = base.Toolbox()

        hof = tools.HallOfFame(hofn)
        self.hof = hof
        
        stats = tools.Statistics(key=lambda ind: 1.0/ind.fitness.values[0])
        stats.register("avg", np.mean)
        stats.register("std", np.std)
        stats.register("min", np.min)
        stats.register("max", np.max)
        self.stats = stats

        # Our fitness already takes into account all the molecules simultaneously.
        # Therefore, there is no need for a multi-objective optimization.
        creator.create("FitnessMin", base.Fitness, weights=(1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMin)
        toolbox.register("randomize", initFunc)
        toolbox.register("individual", tools.initIterate, creator.Individual,
                toolbox.randomize)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        toolbox.register("mate", cxFunc)
        toolbox.register("mutate", mutFunc)
        toolbox.register("select", selFunc)
        toolbox.register("evaluate", evalFunc)

        self.toolbox = toolbox

    def run(self, nGens, nprocs=1):
        # Start processes
        if (nprocs > 1):
            pool = multiprocessing.Pool(processes=nprocs)
            self.toolbox.register("map", pool.map)
        # Run the GA and store the Logbook
        self.population = self.toolbox.population(self.popSize)
        self.output = algorithms.eaSimple(self.population, self.toolbox, self.cxpb,
                self.mutpb, nGens, stats=self.stats, halloffame=self.hof, verbose=True)

    def getBest(self):
        # Return best multiprofile
        best = self.hof[0]
        mp = multiProfile()
        mp.setOptimizableParameters(slice(0, len(best), 1), list(best))
        mp.minimizeProfiles()
        return mp
