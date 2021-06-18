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

from    random  import  uniform, gauss, lognormvariate, choice
from    abc     import  ABC
import  numpy   as      np

class RandomizerInterface(ABC):
    def random(self): pass

class LogNormalRandomizer(RandomizerInterface):
    def __init__(self, mean, stdev):
        self.mean = mean
        self.stdev = stdev
    def random(self):
        phi = (self.stdev ** 2 + self.mean ** 2) ** 0.5
        mu = np.log(self.mean ** 2 / phi)
        sigma = (np.log(phi ** 2 / self.mean ** 2)) ** 0.5
        samples = lognormvariate(mu, sigma)
        return samples

class GaussianRandomizer(RandomizerInterface):
    def __init__(self, mean, stdev):
        self.mean = mean
        self.stdev = stdev
    def random(self):
        samples = gauss(self.mean, self.stdev)
        return samples

class UniformRandomizer(RandomizerInterface):
    def __init__(self, low, up):
        self.low = low
        self.up  = up
    def random(self):
        samples = uniform(self.low, self.up)
        return samples

class UniformDimDist(RandomizerInterface):
    def __init__(self, low, up):
        self.low = low
        self.up  = up
    def random(self):
        low_exp = np.log10(self.low)
        up_exp  = np.log10(self.up)
        exp = uniform(low_exp, up_exp)
        return 10 ** exp

class SignReverserDecorator(RandomizerInterface):
    def __init__(self, randomizer, pinv):
        self._randomizer = randomizer
        self.pinv = pinv
    def random(self):
        x = choice(range(1, 101))
        r = self._randomizer.random()
        if x <= self.pinv:
            return r
        else:
            return -1.0 * r

class LimiterDecorator(RandomizerInterface):
    def __init__(self, randomizer, min_, max_):
        self._randomizer = randomizer
        self.min_ = min_
        self.max_ = max_
    def random(self):
        rand = self._randomizer.random()
        while (rand < self.min_) and (rand > self.max_):
            rand = self._randomizer.random()
        return rand

def RandomizerFactory(typestr, **kwargs):
    if (typestr == 'uniform'):
        return UniformRandomizer(**kwargs)
    if (typestr == 'uniform_dim'):
        return UniformDimDist(**kwargs)
    if (typestr == 'gaussian'):
        return GaussianRandomizer(**kwargs)
    if (typestr == 'lognormal'):
        return LogNormalRandomizer(**kwargs)
