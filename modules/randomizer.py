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

from random import uniform, gauss, lognormvariate, randint, choice
import numpy as np

class parameterRandomizer (object):

    @staticmethod
    def lognormalDist(mean, stdev, n=1):
        phi = (stdev ** 2 + mean ** 2) ** 0.5
        mu = np.log(mean ** 2 / phi)
        sigma = (np.log(phi ** 2 / mean ** 2)) ** 0.5
        samples = lognormvariate(mu, sigma)
        return samples

    @staticmethod
    def gaussianDist(mean, stdev, n=1):
        samples = gauss(mean, stdev)
        return samples

    @staticmethod
    def uniformDist(low, up, n=1):
        samples = uniform(low, up)
        return samples

    @staticmethod
    def uniformDimDist(low, up, n=1):
        low_exp = np.log10(low)
        up_exp  = np.log10(up)
        exp = uniform(low_exp, up_exp)
        return 10 ** exp

    @staticmethod
    def randomizeDihedral (distflag, min, max, mean, stddev, pinv, dihtype):
        if (dihtype == 'standard'):
            output_phi = 0
            output_k = 0
            output_m = 0
            sample = min - 1
            while (sample < min) or (sample > max):
                if (distflag == 1):
                    sample = parameterRandomizer.uniformDist(min,max,1)
                if (distflag == 2):
                    sample = parameterRandomizer.lognormalDist(mean,stddev,1)
                if (distflag == 3):
                    sample = parameterRandomizer.gaussianDist(mean,stddev,1)
            if (distflag == 4):
                sample = parameterRandomizer.uniformDimDist(min,max,1)
            rand = randint(1,100)
            if rand < pinv:
                output_k = -sample
            else:
                output_k = sample

            output_m = randint(1,6)
            p = randint(0,1)
            output_phi = 0
            if p == 1:
                output_phi = 180.00
            return [output_phi, output_k, output_m]
        elif (dihtype == 'ryckaert'):
            output_k = 0
            sample = min - 1
            if (distflag == 4):
                sample = parameterRandomizer.uniformDimDist(min,max,1)
            else:
                while (sample < min) or (sample > max):
                    if (distflag == 1):
                        sample = parameterRandomizer.uniformDist(min,max,1)
                    if (distflag == 2):
                        sample = parameterRandomizer.lognormalDist(mean,stddev,1)
                    if (distflag == 3):
                        sample = parameterRandomizer.gaussianDist(mean,stddev,1)
            rand = randint(1,100)
            if rand < pinv:
                output_k = -sample
            else:
                output_k = sample
            return output_k
            
    @staticmethod
    def randomizeDihedralMslots(distflag, min, max, mean, stddev, pinv, dihtype):
        # Obtain random parameters as in the old versions.
        parameters = parameterRandomizer.randomizeDihedral(distflag, min, max, mean, stddev, pinv, dihtype)
        return parameters

    @staticmethod
    def randomizeLJ (distflag_c6, min_c6, max_c6, mean_c6, stddev_c6,
           distflag_c12, min_c12, max_c12, mean_c12, stddev_c12):
        # c6
        if (distflag_c6 == 4):
            c6 = parameterRandomizer.uniformDimDist(min_c6,max_c6,1)
        else:
            c6 = min_c6 - 1
            while (c6 < min_c6) or (c6 > max_c6):
                if (distflag_c6 == 1):
                    c6 = parameterRandomizer.uniformDist(min_c6, max_c6, 1)
                if (distflag_c6 == 2):
                    c6 = parameterRandomizer.lognormalDist(mean_c6, stddev_c6, 1)
                if (distflag_c6 == 3):
                    c6 = parameterRandomizer.gaussianDist(mean_c6, stddev_c6, 1)
        # c12
        if (distflag_c12 == 4):
            c12 = parameterRandomizer.uniformDimDist(min_c12,max_c12,1)
        else:
            c12 = min_c12 - 1
            while (c12 < min_c12) or (c12 > max_c12):
                if (distflag_c12 == 1):
                    c12 = parameterRandomizer.uniformDist(min_c12, max_c12, 1)
                if (distflag_c12 == 2):
                    c12 = parameterRandomizer.lognormalDist(mean_c12, stddev_c12, 1)
                if (distflag_c12 == 3):
                    c12 = parameterRandomizer.gaussianDist(mean_c12, stddev_c12, 1)
        return [c6, c12]
