from    .energy_force_vectorized  import  MMCalculator
from    .minim                    import  steepestDescentsMinimizer,  conjugateGradientMinimizer
from    .configuration            import  ensemble
from    .randomizer               import  parameterRandomizer
from    .opts                     import  *
from    random                    import  choice,random
import  numpy                     as      np

class profile (object):
    
    def __init__ (self, stpData, trajFile,
                  scanFirst, scanStep, scanLast, restrConst,
                  emAlgo, emDX0, emDXM, emDele, emSteps):
            self.enerProfile        = []
            self.ensemble           = ensemble([])
            self.mmCalc             = MMCalculator()
            if stpData['defaults']['comb-rule'] == 2:
                self.mixType = 'arithmetic'
            else:
                self.mixType = 'geometric'
            self.resetMinim(emAlgo, emDX0, emDXM, emDele, emSteps)
            self.refDih             = stpData['propers'][0][stpData['refdihedral']][:4] # (ai, aj, ak, al) quadruple
            self.restrConst         = restrConst
            self.refPhi             = [scanFirst + i*scanStep for i in range(int((scanLast-scanFirst)/scanStep)+1)]
            self.optDihIdxs         = stpData['optdihedrals'].copy()
            self.optType            = stpData['opttype']
            if self.optType == 'atom':
                self.optAtomIdxs = stpData['optatoms'].copy()
            elif self.optType == 'pair':
                self.optPairIdxs = stpData['optpairs'].copy() 
            # initialize data
            self.ensemble.readFromTrajectory(trajFile)
            self.mmCalc.createFromStpDictionary(stpData)

    def resetMinim(self, emAlgo, emDX0, emDXM, emDele, emSteps):
        if (emAlgo == 1): # SDEM
            self.emAlgo = steepestDescentsMinimizer(dx0=emDX0, maxSteps=emSteps, dxm=emDXM, dele=emDele, mmCalc=self.mmCalc)
        elif (emAlgo == 2): # CGEM
            self.emAlgo = conjugateGradientMinimizer(dx0=emDX0, nsteps=emSteps, dxm=emDXM, prec=emDele, calc=self.mmCalc)

    def setDihedralParameters (self, whichTorsion, phi=None, k=None):
        if (optOpts.dihType == 'standard'):
            if (self.mmCalc.dihedralTerms.getType() == 'standard'):
                self.mmCalc.setOptDihedralParameters(optOpts.optTors[whichTorsion]+1, phi, k)
            else:
                    raise Exception(".inp requests standard dihedrals, but .stp specifies another type")
        elif (optOpts.dihType == 'ryckaert'):
            for idx in (self.optDihIdxs):
                if (self.mmCalc.dihedralTerms.getType() == 'ryckaert'):
                    self.mmCalc.setDihedralParametersRyck(idx, optOpts.optTors[whichTorsion], k)
                else:
                    raise Exception(".inp requests Ryckaert-Bellemanns dihedrals, but .stp specifies another type")
        else:
            raise Exception()

    def setLJParameters (self, cs6=None, cs12=None):
        # If an input parameter is None, it will keep its current value.
        # This behavior is implemented in the setLJParametersForAtom methods.
        if (self.optType == 'atom'):
            for idx in self.optAtomIdxs:
                self.mmCalc.setLJParametersForAtom(idx, cs6, cs12, self.mixType)
        elif (self.optType == 'pair'):
            for idx in self.optPairIdxs:
                self.mmCalc.setLJParametersForPair(idx, cs6, cs12)
        else:
            raise ValueError("setLJParameters for {} is not supported.".format(self.optType))

    def minimizeProfile (self):
        self.enerProfile = []
        for i,phi in enumerate(self.refPhi):
            # set restraints in a push-pop fashion


            self.mmCalc.pushDihedralRestraint (*self.refDih, phi, self.restrConst)
            # minimize
            self.emAlgo.run(self.ensemble[i])
            self.enerProfile.append(self.emAlgo.getUnrestrainedEnergy())

            self.mmCalc.popDihedralRestraint()
        # shift to zero
        self.enerProfile = np.array(self.enerProfile) - np.min(self.enerProfile)
        return True

    def rmsdToData (self, data, wei):
        rmsd = 0
        wei_sum = 0
        n = len(self.enerProfile)
        for i in range(n):
            rmsd += wei[i] * np.abs(self.enerProfile[i] - data[i])
            wei_sum += wei[i]
        rmsd /= wei_sum
        return rmsd

    def saveProfile (self, fn):
        np.savetxt(fn, self.enerProfile)

    def saveTraj (self, fn):
        self.ensemble.writeToFile(fn)

class multiProfile (object):

    def __init__ (self):

        self.profiles = []
        self.cs6      = 0.0
        self.cs12     = 0.0
        self.phi      = [0.0] * len(optOpts.optTors)
        self.k        = [0.0] * len(optOpts.optTors)
        self.m        = optOpts.optTors + 1

        for i in range(optOpts.nSystems):
            trajFile = cmdlineOpts.trajFiles[i]
            stpData  = optOpts.stpData[i]
            if (optOpts.dihType == 'standard'):
                self.profiles.append(profile(stpData, trajFile,
                    dihrestrOpts.start, dihrestrOpts.step, dihrestrOpts.last, dihrestrOpts.k,
                    minimOpts.minimType, minimOpts.dx0, minimOpts.dxm, minimOpts.dele, minimOpts.maxSteps))
            if (optOpts.dihType == 'ryckaert'):
                self.profiles.append(profile(stpData, trajFile,
                    dihrestrOpts.start, dihrestrOpts.step, dihrestrOpts.last, dihrestrOpts.k,
                    minimOpts.minimType, minimOpts.dx0, minimOpts.dxm, minimOpts.dele, minimOpts.maxSteps))

    def __getitem__ (self, i):
        return self.profiles[i]

    def resetMinim(self, emAlgo, emDX0, emDXM, emDele, emSteps):
        for profile in self.profiles:
            profile.resetMinim(emAlgo, emDX0, emDXM, emDele, emSteps)

    def setDihedralParameters (self, whichTorsion, phi=None, k=None):
        if (optOpts.dihType == 'standard'):
            if phi is not None:
                self.phi[whichTorsion] = phi
            if k is not None:
                self.k[whichTorsion] = k
            for profile in self.profiles:
                profile.setDihedralParameters (whichTorsion, phi, k)
        if (optOpts.dihType == 'ryckaert'):
            if k is not None:
                self.k[whichTorsion] = k
            for profile in self.profiles:
                profile.setDihedralParameters (whichTorsion, phi, k)

    def setLJParameters (self, cs6=None, cs12=None):
        if cs6 is not None:
            self.cs6 = cs6
        if cs12 is not None:
            self.cs12 = cs12
        if (optOpts.nLJ == 0):
            self.cs6 = 0
            self.cs12 = 0
            cs6 = None
            cs12 = None
        elif (optOpts.nLJ == -1):
            cs12 = None
            self.cs12 = 0
        elif (optOpts.nLJ == -2):
            cs6 = None
            self.cs6 = 0
        for profile in self.profiles:
            profile.setLJParameters (cs6, cs12)

    def minimizeProfiles (self):
        for profile in self.profiles:
            if not (profile.minimizeProfile()):
                return False
        return True

    def genRandDihedralParameters (self, m):
        parameters = parameterRandomizer.randomizeDihedralMslots(
            randOpts.torsDist, randOpts.torsMin, randOpts.torsMax, randOpts.torsMean, randOpts.torsStddev,
            randOpts.torsPinv, optOpts.dihType)
        if (optOpts.dihType == 'standard'):
            if (m is not None):
                # force 'm' and 'phi' values
                parameters[0] = randOpts.mslots_phi
                parameters[2] = m
        elif (optOpts.dihType == 'ryckaert'):
            pass
        return parameters

    def genRandLJParameters (self):
        return parameterRandomizer.randomizeLJ(
            randOpts.cs6Dist, randOpts.cs6Min, randOpts.cs6Max, randOpts.cs6Mean, randOpts.cs6Stddev,
            randOpts.cs12Dist, randOpts.cs12Min, randOpts.cs12Max, randOpts.cs12Mean, randOpts.cs12Stddev)

    def randomizeDihedralParameters (self):
        for i,k in enumerate(optOpts.optTors):
            randomPars = self.genRandDihedralParameters(k+1)
            if optOpts.dihType == 'standard':
                self.setDihedralParameters(i, k=randomPars[1], phi=randomPars[0])
            if optOpts.dihType == 'ryckaert':
                self.setDihedralParameters(i, k=randomPars)
                
    def randomizeLJParameters (self):
        if (optOpts.nLJ != 0):
            randomPars = self.genRandLJParameters()
            if (optOpts.nLJ == -1):
                randomPars[1] = None
            if (optOpts.nLJ == -2):
                randomPars[0] = None
            self.setLJParameters (*randomPars)

    def mutateLJ (self):
        if (optOpts.nLJ != 0):
            mutPars = self.genRandLJParameters()
            if (optOpts.nLJ == -1):
                chosen = 0
            if (optOpts.nLJ == -2):
                chosen = 1
            if (optOpts.nLJ == 1):
                chosen  = choice(range(len(mutPars)))
            mutPars = [x if i == chosen else None for i,x in enumerate(mutPars)]
            self.setLJParameters (*mutPars)

    def mutateDihedral (self):
        if (optOpts.nTors == 0):
            return
        chosenTorsion = choice(optOpts.optTors)
        mutPars       = self.genRandDihedralParameters(chosenTorsion+1)
        if (optOpts.dihType == 'standard'):
            # i == 1 means you only get 'k'
            mutPars = [x if i == 1 else None for i,x in enumerate(mutPars)]
            self.setDihedralParameters (np.argwhere(optOpts.optTors == chosenTorsion)[0][0], k=mutPars[1], phi=mutPars[0])
        if (optOpts.dihType == 'ryckaert'):
            self.setDihedralParameters (np.argwhere(optOpts.optTors == chosenTorsion)[0][0], k=mutPars)

    @staticmethod
    def uniformCrossover (father, mother, childOne, childTwo):
        chosen   =  choice([0,1])
        if (chosen == 0):
            cs6one   =  father.cs6
            cs6two   =  mother.cs6
        else:
            cs6one   =  mother.cs6
            cs6two   =  father.cs6
        chosen  =  choice([0,1])
        if (chosen == 0):
            cs12one  =  father.cs12
            cs12two  =  mother.cs12
        else:
            cs12one  =  mother.cs12
            cs12two  =  father.cs12
        childOne.setLJParameters (cs6one, cs12one)
        childTwo.setLJParameters (cs6two, cs12two)
        for i in range(len(optOpts.optTors)):
            chosen  =  choice([0,1])
            if (chosen == 0):
                k_1    =  father.k[i]
                k_2    =  mother.k[i]
            else:
                k_1    =  mother.k[i]
                k_2    =  father.k[i]
            chosen  =  choice([0,1])
            if (chosen == 0):
                m_1    =  father.m[i]
                m_2    =  mother.m[i]
            else:
                m_1    =  mother.m[i]
                m_2    =  father.m[i]
            chosen  =  choice([0,1])
            if (chosen == 0):
                phi_1  =  father.phi[i]
                phi_2  =  mother.phi[i]
            else:
                phi_1  =  mother.phi[i]
                phi_2  =  father.phi[i]

            if (optOpts.dihType == 'standard'):
                childOne.setDihedralParameters(i, phi_1, k_1)
                childTwo.setDihedralParameters(i, phi_2, k_2)
            elif (optOpts.dihType == 'ryckaert'):
                childOne.setDihedralParameters(i, k=k_1)
                childTwo.setDihedralParameters(i, k=k_2)
        childOne.minimizeProfiles()
        childTwo.minimizeProfiles()
        return [childOne, childTwo]

    @staticmethod
    def averageCrossover (father, mother, childOne, childTwo):
        cs6one   =  0.50     *  (father.cs6   +  mother.cs6)
        cs12one  =  0.50     *  (father.cs12  +  mother.cs12)
        childOne.setLJParameters (cs6one, cs12one)
        childTwo.setLJParameters (cs6one, cs12one)
        for i in range(len(optOpts.optTors)):
            k_1    =  0.50 * (father.k[i] + mother.k[i])
            k_2    =  0.50 * (father.k[i] + mother.k[i])
            chosen  =  choice([0,1])
            if (chosen == 0):
                m_1    =  father.m[i]
                m_2    =  mother.m[i]
            else:
                m_1    =  mother.m[i]
                m_2    =  father.m[i]
            chosen  =  choice([0,1])
            if (chosen == 0):
                phi_1  =  father.phi[i]
                phi_2  =  mother.phi[i]
            else:
                phi_1  =  mother.phi[i]
                phi_2  =  father.phi[i]
            if (optOpts.dihType == 'standard'):
                childOne.setDihedralParameters(i, phi_1, k_1)
                childTwo.setDihedralParameters(i, phi_2, k_2)
            elif (optOpts.dihType == 'ryckaert'):
                childOne.setDihedralParameters(i, k=k_1)
                childTwo.setDihedralParameters(i, k=k_2)
        childOne.minimizeProfiles()
        childTwo.minimizeProfiles()
        return [childOne, childTwo]

    @staticmethod
    def arithmeticCrossover (father, mother, childOne, childTwo):
        r = random()
        cs6one   = r * father.cs6  + (1 - r) *  mother.cs6
        cs12one  = r * father.cs12 + (1 - r) *  mother.cs12
        cs6two   = (1 - r) * father.cs6  + r * mother.cs6
        cs12two  = (1 - r) * father.cs12 + r * mother.cs12
        childOne.setLJParameters (cs6one, cs12one)
        childTwo.setLJParameters (cs6two, cs12two)
        for i in range(len(optOpts.optTors)):
            k_1    =  r * father.k[i] + (1 - r) * mother.k[i]
            k_2    =  r * mother.k[i] + (1 - r) * father.k[i]
            chosen  =  choice([0,1])
            if (chosen == 0):
                m_1    =  father.m[i]
                m_2    =  mother.m[i]
            else:
                m_1    =  mother.m[i]
                m_2    =  father.m[i]
            chosen  =  choice([0,1])
            if (chosen == 0):
                phi_1  =  father.phi[i]
                phi_2  =  mother.phi[i]
            else:
                phi_1  =  mother.phi[i]
                phi_2  =  father.phi[i]
            if (optOpts.dihType == 'standard'):
                childOne.setDihedralParameters(i, phi_1, k_1)
                childTwo.setDihedralParameters(i, phi_2, k_2)
            elif (optOpts.dihType == 'ryckaert'):
                childOne.setDihedralParameters(i, k=k_1)
                childTwo.setDihedralParameters(i, k=k_2)
        childOne.minimizeProfiles()
        childTwo.minimizeProfiles()
        return [childOne, childTwo]

    @staticmethod
    def heuristicCrossover (bestParent, worstParent, childOne, childTwo, limit=500):
        r = random()

        nlimits = 0
        while True:
            nlimits += 1
            cs6_1 = bestParent.cs6 + r * (bestParent.cs6 - worstParent.cs6)
            cs12_1 = bestParent.cs12 + r * (bestParent.cs12 - worstParent.cs12)
            cs6_2  = r * bestParent.cs6  + (1 - r) *  worstParent.cs6
            cs12_2 = r * bestParent.cs12 + (1 - r) *  worstParent.cs12
            if (cs6_1 > 0) and (cs6_2 > 0) and (cs12_1 > 0) and (cs12_2 > 0):
                break
            else:
                r *= 0.50
            if nlimits >= limit:
                raise Exception("Heuristic crossover failed with cs6 = {:e} and cs12 = {:e}.".format(cs6_1, cs12_1))
        childOne.setLJParameters (cs6_1, cs12_1)
        childTwo.setLJParameters (cs6_2, cs12_2)
        for i in range(len(optOpts.optTors)):
            k_1 = bestParent.k[i] + r * (bestParent.k[i] - worstParent.k[i])
            k_2 = bestParent.k[i] * r + (1 - r) * worstParent.k[i]
            chosen  =  choice([0,1])
            if (chosen == 0):
                m_1    =  bestParent.m[i]
                m_2    =  worstParent.m[i]
            else:
                m_1    =  worstParent.m[i]
                m_2    =  bestParent.m[i]
            chosen  =  choice([0,1])
            if (chosen == 0):
                phi_1  =  bestParent.phi[i]
                phi_2  =  worstParent.phi[i]
            else:
                phi_1  =  worstParent.phi[i]
                phi_2  =  bestParent.phi[i]
            if (optOpts.dihType == 'standard'):
                childOne.setDihedralParameters(i, phi_1, k_1)
                childTwo.setDihedralParameters(i, phi_2, k_2)
            elif (optOpts.dihType == 'ryckaert'):
                childOne.setDihedralParameters(i, k=k_1)
                childTwo.setDihedralParameters(i, k=k_2)
        childOne.minimizeProfiles()
        childTwo.minimizeProfiles()
        return [childOne, childTwo]

    def randomize (self):
        self.randomizeDihedralParameters()
        self.randomizeLJParameters()
        return self.minimizeProfiles()

    def mutate (self):
        if (optOpts.nTors == 0):
            self.mutateLJ()
        elif (optOpts.nLJ == 0):
            self.mutateDihedral()
        else:
            if (choice([0,1]) == 0):
                self.mutateDihedral()
            else:
                self.mutateLJ()
        self.minimizeProfiles()
        return self

    def rmsdToData (self):
        rmsd = 0
        wei_sum = 0
        for i,profile in enumerate(self.profiles):
            n = len(profile.enerProfile)
            for j in range(n):
                rmsd += optOpts.weiData[i][j] * ((profile.enerProfile[j] - optOpts.refData[i][j]) ** 2)
                wei_sum += optOpts.weiData[i][j]
        rmsd /= wei_sum
        rmsd = np.sqrt(rmsd)
        return rmsd

    def saveProfile (self, prefix):
        for i,profile in enumerate(self.profiles):
            fn = prefix + "_{0}.dat".format(i+1)
            profile.saveProfile(fn)

    def saveTraj (self, prefix, ext):
        for i,profile in enumerate(self.profiles):
            fn = prefix + "_{0}.{1}".format(i+1, ext)
            profile.saveTraj(fn)

    def saveParameters (self, fn):
        with open(fn,'w') as fp:
            if (optOpts.nLJ != 0):
                if (optOpts.nLJ == -2):
                    fp.write("# {:>16}\n".format("CS12"))
                    fp.write("{:>18.7e}\n".format(self.cs12)) 
                if (optOpts.nLJ == -1):
                    fp.write("# {:>16}\n".format("CS6"))
                    fp.write("{:>18.7e}\n".format(self.cs6)) 
                if (optOpts.nLJ == 1):
                    fp.write("# {:>16}{:>18}\n".format("CS6", "CS12"))
                    fp.write("{:>18.7e}{:>18.7e}\n".format(self.cs6, self.cs12))
                    
            if (optOpts.nTors != 0):
                if (optOpts.dihType == 'standard'):
                    fp.write("# {:>16}{:>18}{:>5}\n".format("PHI", "K", "M"))
                    j = 0
                    for i in range(6):
                        if i in (optOpts.optTors):
                            fp.write("{:>18.2f}{:>18.7f}{:>5}\n".format(self.phi[j], self.k[j], self.m[j]))
                            j += 1
                        else:
                            fp.write("{:>18.2f}{:>18.7f}{:>5}\n".format(0, 0, i+1))
                elif (optOpts.dihType == 'ryckaert'):
                    fp.write("#{:>17}{:>18}{:>18}{:>18}{:>18}{:>18}\n".format('C0','C1','C2','C3','C4','C5'))
                    for i in range(6):
                        if (i in optOpts.optTors):
                            fp.write("{:>18.7e}".format(self.k[np.argwhere(optOpts.optTors == i)[0][0]]))
                        else:
                            fp.write("{:>18.7e}".format(0.0))
                    fp.write("\n")
                    
            fp.write("#{:>17}\n".format("WRMSD"))
            fp.write("{:>18.7e}\n".format(self.rmsdToData()))

