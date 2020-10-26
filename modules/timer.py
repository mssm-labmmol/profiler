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

