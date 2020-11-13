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

# -*- coding: utf-8 -*-
import copy
from .opts import randOpts
from random import uniform, random, randint
from numpy import min,max, exp
from multiprocessing import Pool
from sys import stdout

class VBGA:
    def __init__(self, individualClass, fitnessMethod, 
                       selectionMethod, numSelect,     numElitized,
                       crossoverMethod, crossoverRate,
                       mutationMethod,  mutationRate,
                       populationSize,
                       numProc,
                       writer):
        self.populationSize = populationSize
        self.individualClass = individualClass
        self.fitnessMethod = fitnessMethod
        self.selectionMethod = selectionMethod
        self.numSelect = numSelect
        self.numElitized     = numElitized
        self.crossoverMethod = crossoverMethod
        self.crossoverRate = crossoverRate
        self.mutationRate = mutationRate
        self.mutationMethod = mutationMethod
        self.numProc        = numProc
        self.population = []
        self.writer = writer

    # i is simply a dummy variable
    def parallelCreateRandomIndividual(self, i):
        ind = self.individualClass()
        ind.randomize()
        return ind

    def createRandomIndividual(self):
        ind = self.individualClass()
        ind.randomize()
        return ind

    def createInitialPopulation(self):
        print("Creating initial population. This may take a while...", file=stdout)
        # only use multiprocessing if necessary, otherwise there might be an overhead
        if (self.numProc > 1):
            with Pool(self.numProc) as p:
                self.population = p.map(self.parallelCreateRandomIndividual,
                                        range(self.populationSize))
        else:
            self.population = [self.createRandomIndividual() for i in range(self.populationSize)]
        print ("Population complete.", file=stdout)
    
    def evaluateIndividual(self, individual):
        individual.fitValue = self.fitnessMethod(individual)

    def getFitValue(self, individual):
        return individual.fitValue

    def evaluatePopulation(self):
        for individual in self.population:
            self.evaluateIndividual(individual)

    def normalizationFitness(self):
        sumFitness = 0
        for individual in self.population:
            sumFitness += individual.fitValue
        for individual in self.population:
            individual.fitValue = individual.fitValue / sumFitness

    def inverseNormalizationFitness(self):
        sumFitness = 0
        for individual in self.population:
            sumFitness += 1.0 / individual.fitValue
        for individual in self.population:
            individual.fitValue = 1.0 / (individual.fitValue * sumFitness)

    def hardElitism(self):
        """
        Puts the NTEL best individuals in an array.
        """
        self.population.sort(key=self.getFitValue)
        elitized = self.population[:self.numElitized]
        return elitized

    def applySelectionMethod(self):
        if (self.numSelect % 2 == 0):
            selectedPopulation = self.selectionMethod(self.population, selectionSize=self.numSelect)
        else:
            selectedPopulation = self.selectionMethod(self.population, selectionSize=self.numSelect+1)
        return selectedPopulation

    def averageFitness(self):
        sum = 0
        for individual in self.population:
            sum += individual.fitValue
        sum /= len(self.population)
        return sum

    def bestFitness(self):
        max = 100000
        for individual in self.population:
            if individual.fitValue < max:
                max = individual.fitValue
        return max

    def parallelCrossoverAndMutation(self, parents):
        father = parents[0]
        mother = parents[1]
        childOne = copy.deepcopy(mother)
        childTwo = copy.deepcopy(father)
        dice = uniform(0, 100)
        dice_mut1 = uniform(0, 100)
        dice_mut2 = uniform(0, 100)
        if dice < self.crossoverRate:
            if (father.fitValue >= mother.fitValue):
                children = self.crossoverMethod(father, mother, childOne, childTwo)
            else:
                children = self.crossoverMethod(mother, father, childOne, childTwo)
            if dice_mut1 < self.crossoverRate:
                self.mutationMethod(children[0])
            if dice_mut2 < self.crossoverRate:
                self.mutationMethod(children[1])
            return children
        else:
            if dice_mut1 < self.crossoverRate:
                self.mutationMethod(childOne)
            if dice_mut2 < self.crossoverRate:
                self.mutationMethod(childTwo)
            return (childOne, childTwo)

    def applyCrossoverAndMutationMethod(self, selected):
        populationCross = []
        parents = [[selected[i], selected[i+1]] for i in range(0, len(selected), 2)]
        # only use multiprocessing if necessary, otherwise there might be an overhead
        if (self.numProc > 1):       
            with Pool(self.numProc) as p:
                childs = p.map(self.parallelCrossoverAndMutation, parents)
        else:
            childs = map(self.parallelCrossoverAndMutation, parents)
        # put childs in population cross
        for brothers in childs:
            populationCross.append(brothers[0])
            populationCross.append(brothers[1])
        return populationCross[:self.numSelect]

    def getBest(self):
        min = 100000
        best = None
        for individual in self.population:
            if individual.fitValue < min:
                best = individual
                min = individual.fitValue
        return best

    def run(self, maxIterations=100):
        self.createInitialPopulation()
        self.evaluatePopulation()
        savePreviousFitness = 100000
        bestFitness = self.bestFitness()
        avFitness = self.averageFitness()
        self.writer.write(-1, self.population, avFitness, bestFitness)
        for i in range(maxIterations):
            print("Running generation {0}... ".format(i+1), file=stdout, end='')
            elitizedInds = self.hardElitism()
            self.inverseNormalizationFitness()
            selected = self.applySelectionMethod()
            populationCross = self.applyCrossoverAndMutationMethod(selected)
            self.population = elitizedInds + populationCross
            print("Done.",file=stdout)
            self.evaluatePopulation()
            bestFitness = self.bestFitness()
            avFitness = self.averageFitness()
            self.population.sort(key=self.getFitValue)
            self.writer.write(i, self.population, avFitness, bestFitness)
            if(bestFitness<savePreviousFitness):
                print("New best weighted RMSD => {0}".format(bestFitness))
                savePreviousFitness = bestFitness

# This one is bad.
def DefaultSelectionMethod (population, selectionSize=10):
    population.sort(key=lambda x: x.fitValue, reverse=True)
    return population[:selectionSize]  # return 0,1,2....selectionSize

def RoulettSelectionMethod(population, selectionSize=10):  # FitValue must be normalized
    select = []
    while len(select) < selectionSize:
        roulettPoint = random()  # (0.0, 1.0)
        fitPoint = 0
        for individual in population:
            fitPoint +=  individual.fitValue 
            if roulettPoint < fitPoint:
                select.append(individual)
                break
    return select  # return 0,1,2....selectionSize

def TournamentSelectionMethod(population, tournamentSize=2, selectionSize=10):
    select = []
    while len(select) < selectionSize:
        bestFitness = -1
        bestInd = None
        for i in range(tournamentSize):
            ind = population[randint(0, len(population) - 1)]
            if (ind.fitValue > bestFitness):
                bestFitness = ind.fitValue
                bestInd = ind
        select.append(bestInd)
    return select

def RankSelectionMethod(population, selectionSize=10):
    select = []
    popSize = len(population)
    normFactor = PartialHarmonicSum(1, popSize)
    while len(select) < selectionSize:
        roulettPoint = random()  # (0.0, 1.0)
        fitPoint = 0
        for i in range(popSize):
            fitPoint += 1.0 / ((i + 1) * normFactor)
            if roulettPoint < fitPoint:
                select.append(population[i])
                break
    return select

def PartialHarmonicSum(low, up):
    harmSum = 0
    for i in range(low, up + 1):
        harmSum += i
    return harmSum
