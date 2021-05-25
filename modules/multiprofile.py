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

from .energy_force_vectorized import MMCalculator
from    .minim                    import  steepestDescentsMinimizer,  conjugateGradientMinimizer
from    .configuration            import  ensemble
from .opts import cmdlineOpts, optOpts, minimOpts, dihrestrOpts
import  numpy                     as      np
from abc import ABC
from .readopts import MaskLooper, IndexListCreator

class MaybeConditionedParameter:
    def __init__(self, value, cond_func=None):
        self.value = value
        self.cond_func = cond_func
    def check(self):
        if self.cond_func is None:
            return True
        return self.cond_func(self.value)
    def set(self, value):
        self.value = value
    def get_string(self):
        return self.id_string
    def __float__(self):
        return float(self.value)
    def __repr__(self):
        return "{:>8} = {:<18.7e}\n".format(self.get_string(), self.value)
                 
class NonNegativeConditionedParameter(MaybeConditionedParameter):
    @staticmethod
    def positive_check(value):
        return value >= 0
    def __init__(self, value):
        super().__init__(value, cond_func=NonNegativeConditionedParameter.positive_check)
        self.id_string = "nnpar"

class DihedralForceConstant(MaybeConditionedParameter):
    def __init__(self, value, m=None):
        super().__init__(value)
        if (m is not None):
            self.m = m
        self.id_string = "k"
    def __repr__(self):
        return "{:>8}_{} = {:<18.7e}\n".format(self.get_string(), self.m, self.value)

class DihedralPhase(MaybeConditionedParameter):
    def __init__(self, value, m=None):
        super().__init__(value)
        if (m is not None):
            self.m = m
        self.id_string = "phi"
    def __repr__(self):
        return "{:>8}_{} = {:<18.7e}\n".format(self.get_string(), self.m, self.value)

class LJC6(NonNegativeConditionedParameter):
    def __init__(self, value):
        super().__init__(value)
        self.id_string = "cs6"

class LJC12(NonNegativeConditionedParameter):
    def __init__(self, value):
        super().__init__(value)
        self.id_string = "cs12"

def ParameterFactory(typestr, m=None):
    if (typestr == 'c6'):
        return LJC6(0.0)
    elif (typestr == 'c12'):
        return LJC12(0.0)
    elif (typestr == 'k'):
        return DihedralForceConstant(0.0, m)
    elif (typestr == 'phi'):
        return DihedralPhase(0.0, m)

class IOptimizableParameters(ABC):

    def get(self, k):
        return float(self._data[k])

    def get_object(self, k):
        return self._data[k]

    def set(self, k, val):
        self._data[k].set(val)

    def set_type(self, type_, vals):
        i = 0
        for t, v in zip(self._types, vals):
            if (t == type_):
                self._data[i].set(v)
            i += 1

    def get_type(self, type_):
        out = []
        for t, v in zip(self._types, self._data):
            if (t == type_):
                out.append(float(v))
        return out

    def get_type_objects(self, type_):
        out = []
        for t, v in zip(self._types, self._data):
            if (t == type_):
                out.append(v)
        return out

    def get_type_as_dict(self, type_):
        out = []
        for t, v in zip(self._types, self._data):
            if (t == type_):
                out.append( {v.get_string(): float(v)} )
        return out

    def to_list(self):
        return [float(x) for x in self._data]

    def is_unphysical(self):
        for x in self._data:
            if not (x.check()):
                return True
        return False

    def writeToStream(self, fp):
        last_type = -1
        for x, t in zip(self._data, self._types):
            if (t != last_type):
                fp.write("Type {}\n".format(t))
                last_type = t
            fp.write(x.__repr__())
         # if (optOpts.nLJ != 0):
         #     if (optOpts.nLJ == -2):
         #         fp.write("# {:>16}\n".format("CS12"))
         #         fp.write("{:>18.7e}\n".format(self.cs12)) 
         #     if (optOpts.nLJ == -1):
         #         fp.write("# {:>16}\n".format("CS6"))
         #         fp.write("{:>18.7e}\n".format(self.cs6)) 
         #     if (optOpts.nLJ == 1):
         #         fp.write("# {:>16}{:>18}\n".format("CS6", "CS12"))
         #         fp.write("{:>18.7e}{:>18.7e}\n".format(self.cs6, self.cs12))
         # if (optOpts.nTors != 0):
         #     if (optOpts.dihType == 'standard'):
         #         fp.write("# {:>16}{:>18}{:>5}\n".format("PHI", "K", "M"))
         #         j = 0
         #         for i in range(6):
         #             if i in (optOpts.optTors):
         #                 fp.write("{:>18.2f}{:>18.7f}{:>5}\n".format(self.phi[j], self.k[j], self.m[j]))
         #                 j += 1
         #             else:
         #                 fp.write("{:>18.2f}{:>18.7f}{:>5}\n".format(0, 0, i+1))
         #     elif (optOpts.dihType == 'ryckaert'):
         #         fp.write("#{:>17}{:>18}{:>18}{:>18}{:>18}{:>18}\n".format('C0','C1','C2','C3','C4','C5'))
         #         for i in range(6):
         #             if (i in optOpts.optTors):
         #                 fp.write("{:>18.7e}".format(self.k[np.argwhere(optOpts.optTors == i)[0][0]]))
         #             else:
         #                 fp.write("{:>18.7e}".format(0.0))
         #         fp.write("\n")
         # fp.write("#{:>17}\n".format("WRMSD"))
         # fp.write("{:>18.7e}\n".format(self.rmsdToData()))

class FullOptimizableParameters(IOptimizableParameters):
    def __init__(self, LJMasks=None, kMasks=None, phiMasks=None):
        self._data, self._types = MaskLooper(ParameterFactory, LJMasks, kMasks, phiMasks)

# End of new content

class profile (object):
    
    def __init__ (self, stpData, trajFile, scanFirst, scanStep, scanLast, restrConst, emAlgo, emDX0, emDXM, emDele, emSteps, emmData=None):
            self.enerProfile        = []
            self.ensemble           = ensemble([])
            self.emmData = emmData
            self.mmCalc             = MMCalculator()
            if stpData['defaults']['comb-rule'] == 2:
                self.mixType = 'arithmetic'
            else:
                self.mixType = 'geometric'
            self.resetMinim(emAlgo, emDX0, emDXM, emDele, emSteps)
            self.refDih             = stpData['propers'][0][stpData['refdihedral']][:4] # (ai, aj, ak, al) quadruple
            self.restrConst         = restrConst
            self.refPhi             = [scanFirst + i*scanStep for i in range(int((scanLast-scanFirst)/scanStep)+1)]
            self.optType            = stpData['opttype']
            # initialize data
            self.ensemble.readFromTrajectory(trajFile)
            self.mmCalc.createFromStpDictionary(stpData)
            # index list creator
            self.indexList = IndexListCreator(self.optType, stpData['optatoms'], stpData['optpairs'], stpData['optdihedrals'])

    def resetMinim(self, emAlgo, emDX0, emDXM, emDele, emSteps):
        if (emAlgo == 1): # SDEM
            self.emAlgo = steepestDescentsMinimizer(dx0=emDX0, maxSteps=emSteps, dxm=emDXM, dele=emDele, mmCalc=self.mmCalc)
        elif (emAlgo == 2): # CGEM
            self.emAlgo = conjugateGradientMinimizer(dx0=emDX0, nsteps=emSteps, dxm=emDXM, prec=emDele, calc=self.mmCalc)

    def resetMMCalcForMinim(self, stpData):
        self.emmData = None
        self.mmCalc = MMCalculator()
        self.mmCalc.createFromStpDictionary(stpData)

    def prepareMinim(self, emAlgo, emDX0, emDXM, emDele, emSteps, stpData):
        self.resetMMCalcForMinim(stpData)
        self.resetMinim(emAlgo, emDX0, emDXM, emDele, emSteps)

    def setDihedralParameters (self, whichOpt, m, phi=None, k=None):
        # FIXME: This needs to be refactored.
        if (optOpts.dihType == 'standard'):
            for idx in range(len(self.indexList.get(whichOpt))): # Note the range(len(
                if (self.mmCalc.dihedralTerms.getType() == 'standard'):
                    self.mmCalc.setOptDihedralParameters(idx, m, phi, k)
                else:
                    raise Exception(".inp requests standard dihedrals, but .stp specifies another type")
        elif (optOpts.dihType == 'ryckaert'):
            for idx in self.indexList.get(whichOpt): # Note the absence of range(len(
                if (self.mmCalc.dihedralTerms.getType() == 'ryckaert'):
                    self.mmCalc.setDihedralParametersRyck(idx, m-1, k)
                else:
                    raise Exception(".inp requests Ryckaert-Bellemanns dihedrals, but .stp specifies another type")
        else:
            raise Exception()

    def setLJParameters (self, type_, cs6=None, cs12=None):
        # If an input parameter is None, it will keep its current value.
        # This behavior is implemented in the setLJParametersForAtom methods.
        indexList = self.indexList.get(type_)
        if (self.optType == 'atom'):
            self.mmCalc.setLJParametersForAtoms(indexList, cs6, cs12, self.mixType)
        elif (self.optType == 'pair'):
            for idx in indexList:
                self.mmCalc.setLJParametersForPair(idx, cs6, cs12)
        else:
            raise ValueError("setLJParameters for {} is not supported.".format(self.optType))

    def setParameters(self, type_, cs6=None, cs12=None, m=None, phi=None, k=None):
        self.setLJParameters(type_, cs6=cs6, cs12=cs12)
        self.setDihedralParameters(type_, m, phi=phi, k=k)

    def minimizeProfile (self, wei=None, veryLargeEnergy=1.0e+05):
        self.enerProfile = []
        for i,phi in enumerate(self.refPhi):
            if (wei is None) or (wei[i] != 0):
                if (self.emmData is None):
                    # set restraints in a push-pop fashion
                    self.mmCalc.pushDihedralRestraint (*self.refDih, phi_0=phi, k=self.restrConst)
                    # minimize
                    self.emAlgo.run(self.ensemble[i])
                    self.mmCalc.popDihedralRestraint()
                    self.enerProfile.append(
                        self.mmCalc.calcForConf(self.ensemble[i], removeRestraintsFromTotal=True, calcForce=False, minim=False)[0]['total']
                    )
                else:
                    raise Exception("In new code versions, this branch shouldn't be executed.")
                    # self.mmCalc.setEmm(self.emmData[i])
                    # self.emAlgo.run(self.ensemble[i])
                    # self.enerProfile.append(self.mmCalc.calcForConf(self.ensemble[i], calcForce=False, minim=False))
            else:
                self.enerProfile.append(veryLargeEnergy)
        # shift to zero mean
        self.enerProfile = np.array(self.enerProfile) - np.mean(self.enerProfile)
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
        #self.m        = [[s + 1 for s in multSwitches] for multSwitches in optOpts.optTors]
        self.optPars  = None  # It will be initialized separately for each dihedral type.

        for i in range(optOpts.nSystems):
            trajFile = cmdlineOpts.trajFiles[i]
            stpData  = optOpts.stpData[i]
            if optOpts.emmData is None:
                emmData = None
            else:
                emmData = optOpts.emmData[i,:]
            if (optOpts.dihType == 'standard'):
                self.profiles.append(profile(stpData, trajFile,
                    dihrestrOpts.start, dihrestrOpts.step, dihrestrOpts.last, dihrestrOpts.k,
                    minimOpts.minimType, minimOpts.dx0, minimOpts.dxm, minimOpts.dele, minimOpts.maxSteps,
                                             emmData=emmData))

                self.optPars = FullOptimizableParameters(LJMasks=optOpts.LJMask, kMasks=optOpts.kMask, phiMasks=optOpts.phiMask)

            if (optOpts.dihType == 'ryckaert'):
                self.profiles.append(profile(stpData, trajFile,
                    dihrestrOpts.start, dihrestrOpts.step, dihrestrOpts.last, dihrestrOpts.k,
                                             minimOpts.minimType, minimOpts.dx0, minimOpts.dxm, minimOpts.dele, minimOpts.maxSteps,
                                             emmData=emmData))

                # ignores mask for optimization of phases
                self.optPars = FullOptimizableParameters(LJMasks=optOpts.LJMask, kMasks=optOpts.kMask)

    def __getitem__ (self, i):
        return self.profiles[i]

    def getOptimizableParameters(self):
        return self.optPars.to_list()

    def areThereUnphysicalParameters(self):
        return self.optPars.is_unphysical()

    def getNumberOfOptimizableParameters(self):
        return len(self.optPars.to_list())

    def setSingleOptimizableParameter(self, k, val):
        self.optPars.set(k, val)
        self.applyParameters()

    def setOptimizableParameters(self, ks, vals):
        if (type(ks) is slice):
            for k, v in zip(range(ks.start or 0, ks.stop or self.getNumberOfOptimizableParameters(), ks.step or 1), vals):
                self.setSingleOptimizableParameter(k, v)
        else:
            self.setSingleOptimizableParameter(ks, vals)

    def prepareMinim(self, emAlgo, emDX0, emDXM, emDele, emSteps):
        for i, profile in enumerate(self.profiles):
            profile.prepareMinim(emAlgo, emDX0, emDXM, emDele, emSteps, optOpts.stpData[i])
        self.applyParameters()

    def applyParameters(self):
        """Applies current optimizable parameters to the member profiles."""
        for profile in self.profiles:
            for i in range(optOpts.nLJ):
                for p in self.optPars.get_type_as_dict(i):
                    profile.setParameters(i, **p)
            for i in range(optOpts.nLJ, optOpts.nLJ + optOpts.nTors):
                for p, obj in zip(self.optPars.get_type_as_dict(i), self.optPars.get_type_objects(i)):
                    profile.setParameters(i, m=obj.m, **p)

    def getNonoptEnergy (self):
        raise NotImplementedError()
        # # Returns a matrix NSYSTEMS x NDIHS with unrestrained energies
        # # containing only the nonoptimized terms. For each system,
        # # data is shifted so that the average value is zero.
        # nonopt_energies = np.zeros((len(self.profiles), len(self.profiles[0].refPhi)))
        # ind_copy = deepcopy(self)
        # # zero out the parameters
        # if (optOpts.nLJ != 0):
        #     ind_copy.setLJParameters(0, 0)
        # for i,k in enumerate(optOpts.optTors):
        #     ind_copy.setDihedralParameters(i, 0.0, 0.0)
        # for i in range(nonopt_energies.shape[0]):
        #     for j in range(nonopt_energies.shape[1]):
        #         nonopt_energies[i,j] = ind_copy.profiles[i].mmCalc.calcForConf(ind_copy.profiles[i].ensemble[j], removeRestraintsFromTotal=True)[0]['total']
        # nonopt_energies[i,:] -= np.mean(nonopt_energies[i,:])
        # return nonopt_energies

    def minimizeProfiles (self, useWei=False):
        if (useWei):
            raise Exception("useWei is dangerous")
            #for i,profile in enumerate(self.profiles):
            #    if not (profile.minimizeProfile(wei=optOpts.weiData[i][:])):
            #        return False
            #return True
        else:
            for i,profile in enumerate(self.profiles):
                if not (profile.minimizeProfile()):
                    return False
            return True

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
            self.optPars.writeToStream(fp)
        fp.write("rmsd = {}".format(self.rmsdToData()))
