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

import numpy as np
import re
from   subprocess import check_output
from   shutil import which
from   io     import StringIO

def raiseError (msg):
    print ("ERROR: " + msg)
    exit(0)

def checkStpExtension (fn):
    if (fn.endswith('.stp')):
        return True
    return False

def createStreamAfterPreprocessing (fn):
    cpp_path = which('cpp')
    if cpp_path is None:
        answer = 'n'
        while not (answer == 'y'):
            print("Error: The stpParser uses the C preprocessor (cpp) to translate")
            print("directives such as #include and #define. However, no path for cpp")
            print("was found. Please, preprocess the file %s manually. If you have already" % fn)
            print("done this, or if your file does not need preprocessing, answer with")
            print("'y'. If you want to quit, send the kill signal Ctrl+C.")
            answer = input("")
        fp = open(fn, 'r')
    else:
        processed_string = check_output([cpp_path, '-P', '-traditional', fn]).decode('utf-8')
        fp = StringIO(processed_string)
    return fp

# goto block '[ stringId ]'
def gotoBlock (stream, stringId):
    for line in stream:
        if re.match(r"^\[ " + re.escape(stringId) + r" ]", line):
            return

# goto next block and return block name
def gotoNextBlock (stream):
    for line in stream:
        m = re.match("\[ (\w+) \]", line)
        if m:
            return m.group(1)
    return -1

# assumes that stream is currently at a '[ bonds ]' block
def readBondsFromStream (stream):
    outputParticles  = []
    outputParameters = []
    for line in stream:
        if (re.match(r"^;", line)):
            continue
        if (len(line.split()) < 2):
            break
        flds = line.split()
        # check for type
        if (int(flds[2])) not in [1,2]:
            raiseError ("all bonds must be type 1 or 2 and given explicitly")
        outputParticles.append( (int(flds[0]), int(flds[1]), int(flds[2])) ) 
        outputParameters.append( [float(flds[3]), float(flds[4])] )
    outputParameters = np.array(outputParameters)
    return (outputParticles, outputParameters)

# assumes that stream is currently at an '[ angles ]' block
def readAnglesFromStream (stream):
    outputParticles  = []
    outputParameters = []
    for line in stream:
        if (re.match(r"^;", line)):
            continue
        if (len(line.split()) < 2):
            break
        flds = line.split()
        # check for type
        if (int(flds[3])) in [1,2]:
            outputParticles.append( (int(flds[0]), int(flds[1]), int(flds[2]), int(flds[3])) ) 
            outputParameters.append( [float(flds[4]), float(flds[5]), 0.0, 0.0] )
        elif int(flds[3]) == 5:
            outputParticles.append( (int(flds[0]), int(flds[1]), int(flds[2]), int(flds[3])) ) 
            outputParameters.append( [float(flds[4]), float(flds[5]), float(flds[6]), float(flds[7])] )
        else:
            raiseError ("all angles must be type 1, 2 or 5 and given explicitly")
    outputParameters = np.array(outputParameters)
    return (outputParticles, outputParameters)

# Assumes that stream is currently at a '[ dihedrals ]' block.
# This function is a bit trickier because the block may contain proper or improper dihedrals.
# To accomodate these two possibilities, the output tuple distinguishes the proper and improper cases.
def readDihedralsFromStream (stream):
    outputParticles  = []
    outputParameters = []
    outputImproperParticles  = []
    outputImproperParameters = []
    outputRestraintParticles = []
    outputRestraintParameters = []
    for line in stream:
        if (re.match(r"^;", line)):
            continue
        if (len(line.split()) < 2):
            break
        line = line.split(';')[0]
        flds = line.split()
        # check for type
        if (int(flds[4]) == 1) or (int(flds[4]) == 9):
            # proper dihedral
            outputParticles.append((int(flds[0]), int(flds[1]), int(flds[2]), int(flds[3]), 1)) # hard-coded 1 to avoid Exception for different dihedral types
            outputParameters.append([float(flds[5]), float(flds[6]), int(flds[7]), 0.0, 0.0, 0.0])
        elif (int(flds[4]) == 3):
            # Ryckaert-Bellemanns
            outputParticles.append((int(flds[0]), int(flds[1]), int(flds[2]), int(flds[3]), int(flds[4]))) 
            outputParameters.append([float(flds[5]), float(flds[6]), float(flds[7]), float(flds[8]), float(flds[9]), float(flds[10])])
        elif (int(flds[4]) == 5):
            # Fourier
            outputParticles.append((int(flds[0]), int(flds[1]), int(flds[2]), int(flds[3]), int(flds[4]))) 
            outputParameters.append([float(flds[5]), float(flds[6]), float(flds[7]), float(flds[8]), 0.0, 0.0])
        elif (int(flds[4]) == 2):
            # improper
            outputImproperParticles.append( (int(flds[0]), int(flds[1]), int(flds[2]), int(flds[3]), int(flds[4])) ) 
            outputImproperParameters.append( [float(flds[5]), float(flds[6])] )
        elif (int(flds[4]) == 4):
            # improper
            outputImproperParticles.append( (int(flds[0]), int(flds[1]), int(flds[2]), int(flds[3]), int(flds[4])) ) 
            outputImproperParameters.append( [float(flds[5]), float(flds[6]), int(flds[7])] )
        elif (int(flds[4]) == -1):
            # dihedral restraint
            outputRestraintParticles.append( (int(flds[0]), int(flds[1]), int(flds[2]), int(flds[3])) ) 
            outputRestraintParameters.append( [float(flds[5]), float(flds[6])] )
        else:
            raiseError ("dihedral must be type 1 (proper), 2 (improper), 3, (Ryckert-Bellemanns), 4 (periodic improper), 5 (Fourier), 9 (multiple propers) or -1 (dihedral restraint)  and given explicitly")
    outputParameters = np.array(outputParameters)
    outputImproperParameters = np.array(outputImproperParameters)
    outputRestraintParameters = np.array(outputRestraintParameters)
    return ((outputParticles, outputParameters), (outputImproperParticles, outputImproperParameters), (outputRestraintParticles, outputRestraintParameters))

# Assumes that stream is currently at a [ defaults ] block.
# Syntax is nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ.
def readDefaultsFromStream(stream):
    outputDict = {}
    allowedDict = {
        'nbfunc': [1],
        'comb-rule': [1, 2, 3],
        'gen-pairs': ['fudge', 'atomic', 'no']}
    for line in stream:
        if (re.match(r"^;", line)):
            continue
        # remove end comment
        line = line.split(';')[0]
        flds = line.split()
        # assign values
        outputDict['nbfunc'] = int(flds.pop(0))
        outputDict['comb-rule'] = int(flds.pop(0))
        outputDict['gen-pairs'] = flds.pop(0)
        outputDict['fudgeLJ'] = float(flds.pop(0))
        outputDict['fudgeQQ'] = float(flds.pop(0))
        break
    # check allowed values
    for flag in ['nbfunc', 'comb-rule', 'gen-pairs']:
        if (outputDict[flag] not in allowedDict[flag]):
            raise Exception("Value {} not allowed for {}.".format(outputDict[flag], flag))
    return outputDict

# Assumes that stream is currently at a [ atoms ] block.
def readAtomsFromStream (stream, comb_rule=1, gen_pairs='atomic'):
    outputParameters = []
    for line in stream:
        if (re.match(r"^;", line)):
            continue
        # remove end comment
        line = line.split(';')[0]
        flds = line.split()
        nflds = len(flds)
        if (nflds < 2):
            break
        if (nflds == 6):
            if (comb_rule == 1):
                outputParameters.append ({'type': flds[0], 
                                          'c6': float(flds[1]), 
                                          'c12': float(flds[2]), 
                                          'cs6': float(flds[3]), 
                                          'cs12': float(flds[4]), 
                                          'q': float(flds[5])})
            elif ((comb_rule == 2) or (comb_rule == 3)):
                # convert sigma-epsilon to c6-c12
                sigma = float(flds[1])
                epsilon = float(flds[2])
                ssigma = float(flds[3])
                sepsilon = float(flds[4])
                outputParameters.append ({'type': flds[0], 
                                          'c6': 4 * epsilon * (sigma**6),
                                          'c12': 4 * epsilon * (sigma ** 12),
                                          'cs6': 4 * sepsilon * (ssigma ** 6),
                                          'cs12': 4 * sepsilon * (ssigma ** 12),
                                          'q': float(flds[5])})
        elif (nflds == 4):
            if (gen_pairs == 'atomic'):
                raise Exception("Can't generate 1-4 pair parameters if 1-4 atomic parameters are not given.")
            else:
                if (comb_rule == 1):
                    outputParameters.append ({'type': flds[0], 
                                              'c6': float(flds[1]), 
                                              'c12': float(flds[2]), 
                                              'cs6': 0.0,
                                              'cs12': 0.0,
                                              'q': float(flds[3])})
                elif ((comb_rule == 2) or (comb_rule == 3)):
                    # convert sigma-epsilon to c6-c12
                    sigma = float(flds[1])
                    epsilon = float(flds[2])
                    outputParameters.append ({'type': flds[0], 
                                              'c6': 4 * epsilon * (sigma**6),
                                              'c12': 4 * epsilon * (sigma ** 12),
                                              'cs6': 0,
                                              'cs12': 0,
                                              'q': float(flds[3])})
    outputParameters = np.array(outputParameters)
    return outputParameters

# Assumes that stream is currently at a special '[ nbpairs ]' block.
def readNonbondedFromStream (stream, atomicParameters, comb_rule=1, gen_pairs='atomic', fudgeLJ=1.0, fudgeQQ=1.0):
    outputParticles  = []
    outputParameters = []
    for line in stream:
        if (re.match(r"^;", line)):
            continue
        line = line.split(';')[0]
        flds = line.split()
        nflds = len(flds)
        if (nflds < 2):
            break
        outputParticles.append( (int(flds[0]), int(flds[1]), int(flds[2])) ) 
        # outputParameters must be calculated from the atomicParameters
        nb_type = int(flds[2]);
        ai_idx  = int(flds[0]) - 1;
        aj_idx  = int(flds[1]) - 1;
        c6      = 0.0
        c12     = 0.0
        qiqj    = 0.0
        qiqj = atomicParameters[ai_idx]['q'] * atomicParameters[aj_idx]['q']
        if (nb_type == 1):
            # this is a standard type, get c6 and c12 from line or atoms
            if (nflds == 3):
                # from atoms
                c6 = np.sqrt(atomicParameters[ai_idx]['c6']) * np.sqrt(atomicParameters[aj_idx]['c6'])
                c12 = np.sqrt(atomicParameters[ai_idx]['c12']) * np.sqrt(atomicParameters[aj_idx]['c12'])
            else:
                # from line
                if (comb_rule == 1):
                    # line contains c6, c12
                    c6 = float(flds[3])
                    c12 = float(flds[4])
                elif ((comb_rule == 2) or (comb_rule == 3)):
                    # line contains sigma, epsilon 
                    sigma = float(flds[3])
                    epsilon = float(flds[4])
                    c6 = 4 * epsilon * (sigma ** 6)
                    c12 = 4 * epsilon * (sigma ** 12)
        elif (nb_type == 2):
            # this is a 1-4 pair
            # charges are fudged
            qiqj *= fudgeQQ
            # Lennard-Jones depends
            if (nflds == 3):
                # read default
                if (gen_pairs == 'atomic'):
                    # mix atomic parameters
                    if (comb_rule == 1) or (comb_rule == 3):
                        c6 = np.sqrt(atomicParameters[ai_idx]['cs6']) * np.sqrt(atomicParameters[aj_idx]['cs6'])
                        c12 = np.sqrt(atomicParameters[ai_idx]['cs12']) * np.sqrt(atomicParameters[aj_idx]['cs12'])
                    elif (comb_rule == 2):
                        c6i = atomicParameters[ai_idx]['cs6']
                        c12i = atomicParameters[ai_idx]['cs12']
                        c6j = atomicParameters[aj_idx]['cs6']
                        c12j = atomicParameters[aj_idx]['cs12']
                        sigmai = (c12i/c6i)**(1.0/6.0)
                        epsiloni = 0.25 * (c6i**2) / c12i
                        sigmaj = (c12j/c6j)**(1.0/6.0)
                        epsilonj = 0.25 * (c6j**2) / c12j
                        sigma = 0.5 * (sigmai + sigmaj)
                        epsilon = np.sqrt(epsiloni * epsilonj)
                        c6 = 4 * epsilon * (sigma ** 6)
                        c12 = 4 * epsilon * (sigma ** 12)
                elif (gen_pairs == 'fudge'):
                    if (comb_rule == 1) or (comb_rule == 3):
                        c6 = fudgeLJ * np.sqrt(atomicParameters[ai_idx]['c6']) * np.sqrt(atomicParameters[aj_idx]['c6'])
                        c12 = fudgeLJ * np.sqrt(atomicParameters[ai_idx]['c12']) * np.sqrt(atomicParameters[aj_idx]['c12'])
                    elif (comb_rule == 2):
                        c6i = atomicParameters[ai_idx]['c6']
                        c12i = atomicParameters[ai_idx]['c12']
                        c6j = atomicParameters[aj_idx]['c6']
                        c12j = atomicParameters[aj_idx]['c12']
                        sigmai = (c12i/c6i)**(1.0/6.0)
                        epsiloni = 0.25 * (c6i**2) / c12i
                        sigmaj = (c12j/c6j)**(1.0/6.0)
                        epsilonj = 0.25 * (c6j**2) / c12j
                        sigma = 0.5 * (sigmai + sigmaj)
                        epsilon = np.sqrt(epsiloni * epsilonj)
                        c6 = fudgeLJ * 4 * epsilon * (sigma ** 6)
                        c12 = fudgeLJ * 4 * epsilon * (sigma ** 12)                        
                elif (gen_pairs == 'no'):
                    raise Exception("Must specify pair parameters if gen-pairs = no.")
            else:
                # from line
                if (comb_rule == 1):
                    # line contains c6, c12
                    c6 = float(flds[3])
                    c12 = float(flds[4])
                elif ((comb_rule == 2) or (comb_rule == 3)):
                    # line contains sigma, epsilon 
                    sigma = float(flds[3])
                    epsilon = float(flds[4])
                    c6 = 4 * epsilon * (sigma ** 6)
                    c12 = 4 * epsilon * (sigma ** 12)
        else:
            # Stop!
            raiseError ("non-bonded pair type must be 1 or 2, not %d" % nb_type)
        outputParameters.append( [c6, c12, qiqj] )
    outputParameters = np.array(outputParameters)
    return (outputParticles, outputParameters)

def readListFromStreamAndSubtractOne (stream):
    out = []
    for line in stream:
        if (re.match(r"^;", line)):
            continue
        if (len(line.split()) < 1):
            break
        out += [int(x)-1 for x in line.split()]
    return out

def readDihedralBlock(stream):
    out = []
    for line in stream:
        if (re.match(r"^;", line)):
            continue
        line = line.split(';')[0]
        flds = line.split()
        if flds == []:
            break
        if len(flds) != 4:
            raise Exception("Wrong specification of dihedral.")
        out.append([int(x) for x in line.split()])
    return out

def readOptpairBlock(stream, nbparticles):
    out = []
    for line in stream:
        if (re.match(r"^;", line)):
            continue
        line = line.split(';')[0]
        flds = line.split()
        if flds == []:
            break
        if len(flds) != 2:
            raise Exception("Wrong specification of optpairs.")
        try:
            out.append(nbparticles.index((int(flds[0]), int(flds[1]), 2)))
        except ValueError:
            raise ValueError("optpair not found.")
    return out    

def parseStpFile (filename, prepareOpt=False):
    defaults = None
    atoms = []
    bondParticles = []
    angleParticles = []
    dihedralParticles = []
    improperParticles = []
    restraintsParticles = []
    nbParticles = []
    bondK = []
    bondL = []
    angleK = []
    angleT = []
    angleKub = []
    angleR13 = []
    dihedralPhi = []
    dihedralK = []
    dihedralM = []
    dihedralC3 = []
    dihedralC4 = []
    dihedralC5 = []
    improperK = []
    improperPhi = []
    improperM   = [] # for periodic
    restraintsK = []
    restraintsPhi = []
    nbC6 = []
    nbC12 = []
    nbQIJ = []
    optAtomsIdxs = []
    optDihIdxs = []
    reflist = []
    optType = None
    optPairIdxs = None
    fp = createStreamAfterPreprocessing(filename)
    while True:
        bn = gotoNextBlock(fp)
        if (bn == -1):
            break
        if (bn == 'defaults'):
            defaults = readDefaultsFromStream(fp)
        if (bn == 'atoms'):
            if (defaults is None):
                raise Exception("[ defaults ] block must come before [ atoms ]")
            atoms = readAtomsFromStream (fp, defaults['comb-rule'], defaults['gen-pairs'])
        elif (bn == "bonds"):
            (particles, pars) = readBondsFromStream (fp)
            bondParticles += particles
            if (pars.shape[0] > 0):
                bondK = np.append(bondK,pars[:,1])
                bondL = np.append(bondL,pars[:,0])
        elif(bn == "angles"):
            (particles,pars) = readAnglesFromStream (fp)
            angleParticles += particles
            if (pars.shape[0] > 0):
                angleK = np.append(angleK,pars[:,1])
                angleT = np.append(angleT,pars[:,0])
                angleKub = np.append(angleKub, pars[:,3])
                angleR13 = np.append(angleR13, pars[:,2])
        elif(bn == "dihedrals"):
            ((particles,pars), (particlesI,parsI), (particlesR, parsR)) = readDihedralsFromStream(fp)
            dihedralParticles += particles
            improperParticles += particlesI
            restraintsParticles += particlesR
            if (pars.shape[0] > 0):
                dihedralPhi       = np.append(dihedralPhi,pars[:,0])
                dihedralK         = np.append(dihedralK,pars[:,1])
                dihedralM         = np.append(dihedralM,pars[:,2])
                dihedralC3         = np.append(dihedralC3,pars[:,3])
                dihedralC4         = np.append(dihedralC4,pars[:,4])
                dihedralC5         = np.append(dihedralC5,pars[:,5])
            if (parsI.shape[0] > 0):
                improperK         = np.append(improperK,parsI[:,1])
                improperPhi       = np.append(improperPhi,parsI[:,0])
                if parsI.shape[1] > 2:
                    improperM = np.append(improperM,parsI[:,2])
                else:
                    improperM = np.append(improperM, [0.0 for i in range(parsI.shape[0])])
            if (parsR.shape[0] > 0):
                restraintsK         = np.append(restraintsK,parsR[:,1])
                restraintsPhi       = np.append(restraintsPhi,parsR[:,0])
        elif(bn == "nbpairs"):
            (particles,pars)  = readNonbondedFromStream (fp, atoms, defaults['comb-rule'], defaults['gen-pairs'], defaults['fudgeLJ'], defaults['fudgeQQ'])
            nbParticles += particles
            if (pars.shape[0] > 0):
                nbC6 = np.append(nbC6,pars[:,0])
                nbC12= np.append(nbC12,pars[:,1])
                nbQIJ = np.append(nbQIJ,pars[:,2])
        elif(bn == "optatoms"):
            if optType is None:
                optAtomsIdxs = readListFromStreamAndSubtractOne (fp)
                optType = 'atom'
                if (defaults['gen-pairs'] == 'no') and (prepareOpt):
                    err_str = "Cannot perform optimization of atomic 1,4 parameters with gen-pairs = no.\n"
                    raise ValueError(err_str)
            else:
                raise ValueError("Can only set one optType.")
        elif(bn == 'optpairs'):
            if optType is None:
                optPairIdxs = readOptpairBlock(fp, nbParticles)
                optType = 'pair'
            else:
                raise ValueError("Can only set one optType.")
        elif(bn == "optdihedrals"):
            optDihIdxs = readDihedralBlock (fp)
            # check if they belong
            dihparticles = [[x[0], x[1], x[2], x[3]] for x in dihedralParticles]
            for dih in optDihIdxs:
                if dih not in dihparticles:
                    raise RuntimeError("Optimized dihedral not found in [ dihedrals ].")
        elif(bn == "refdihedral"):
            reflist = readDihedralBlock(fp)
            # check if it belongs
            dihparticles = [[x[0], x[1], x[2], x[3]] for x in dihedralParticles]
            if reflist[0] not in dihparticles:
                raise RuntimeError("Reference dihedral not found in [ dihedrals ].")
            if (len(reflist) == 0):
                raise RuntimeError("You need to specify one reference dihedral per system!")
            if (len(reflist) > 1):
                raise RuntimeError("You can only list one reference dihedral per system!")
            refDihIdx = reflist[0]
    fp.close()
    if (prepareOpt):
        dihparticles = [[x[0], x[1], x[2], x[3]] for x in dihedralParticles]
        # For each dihedral in [ optdihedrals ], remove multiple instances from stp data.
        for dihAtoms in optDihIdxs:
            idxInPropers = []
            for i,dih in enumerate(dihparticles):
                if (dihAtoms == dih):
                    idxInPropers.append(i)
            if len(idxInPropers) > 1:
                for x in sorted(idxInPropers[1:], reverse=True):
                    del dihedralParticles[x]
                    del dihparticles[x]
                    dihedralPhi = np.delete(dihedralPhi, x)
                    dihedralK = np.delete(dihedralK, x)
                    dihedralM = np.delete(dihedralM, x)
                    dihedralC3 = np.delete(dihedralC3, x)
                    dihedralC4 = np.delete(dihedralC4, x)
                    dihedralC5 = np.delete(dihedralC5, x)
        # Recover indices.
        refDihIdx = dihparticles.index(refDihIdx)
        optDihIdxs = [dihparticles.index(d) for d in optDihIdxs]
    else:
        optAtomsIdxs = []
        optPairIdxs = []
        optDihIdxs = []
        refDihIdx = dihparticles.index(refDihIdx)
    out = {
        'defaults': defaults,
        'atoms': atoms,
        'bonds': (bondParticles, bondK, bondL),
        'angles': (angleParticles, angleK, angleT, angleKub, angleR13),
        'propers': (dihedralParticles, dihedralPhi, dihedralK, dihedralM, dihedralC3, dihedralC4, dihedralC5), 
        'impropers': (improperParticles, improperK, improperPhi, improperM),
        'restraints': (restraintsParticles, restraintsK, restraintsPhi),
        'nb': (nbParticles, nbC6, nbC12, nbQIJ),
        'optatoms': optAtomsIdxs,
        'optpairs': optPairIdxs,
        'optdihedrals': optDihIdxs,
        'refdihedral': refDihIdx, 
        'opttype': optType
    }
    return out
