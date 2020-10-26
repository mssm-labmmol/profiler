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

from time import time
from datetime import timedelta

class Timer (object):

    def __init__ (self, name):

        self.name = name
        self.timeCount = 0
        self.currentCount = 0
        self.isActive  = False
        self.numCalls = 0

    def __del__ (self):
        if (self.isActive):
            print ("Warning: Timer %s is being destroyed, but it is still active." % self.name)
            print ("         It will be automatically deactivated.")
            self.deactivate()
        self.printLog()
        del(self.name)
        del(self.timeCount)
        del(self.currentCount)
        del(self.isActive)

    @staticmethod
    def raiseActiveError ():
        print ("Trying to activate a Timer that is already active.")
        exit()

    @staticmethod
    def raiseInactiveError ():
        print ("Trying to deactivate a Timer that is already inactive.")
        exit()

    def activate (self):
        if (self.isActive):
            Timer.raiseActiveError
        else:
            self.currentCount = time()
            self.isActive = True

    def deactivate (self):
        if not (self.isActive):
            Timer.raiseInactiveError
        else:
            self.timeCount += time() - self.currentCount
            self.numCalls += 1
            self.isActive = False

    def printLog (self):
        print ("Timer %s:" % self.name)
        print ("\tRuntime: %s" % str(timedelta(seconds=self.timeCount)))
        print ("\tCalled %d times." % self.numCalls)

