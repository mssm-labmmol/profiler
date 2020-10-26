#!/usr/bin/env python3

class Option:

    def __init__(self, name, flag, values, default_value):
        self.name  = name
        self.flag  = flag
        self.values = values
        self.value = default_value
        self.default_value = default_value
        self.state = -1

    def advance(self):
        self.state += 1
        if (self.state == len(self.values)):
            self.state = -1
            self.value = self.default_value
            return False
        self.value = self.values[self.state]
        return True

    def putInString(self, string):
        new_string = string.replace(self.flag, str(self.value))
        if (new_string == string):
            raise Exception("Nothing to substitute for flag {}!".format(self.flag))
        return new_string

class FileSweeper:

    def __init__(self, fn, opts):
        self.fn = fn
        self.opts = {}
        for opt in opts:
            self.opts[opt.name] = opt

    def _putCurrentState(self, string):
        new_string = string
        for opt in self.opts.values():
            new_string = opt.putInString(new_string)
        return new_string

    def writeCurrentState(self, fno):
        fi = open(self.fn, 'r')
        string = fi.read()
        fi.close()
        new_string = self._putCurrentState(string)
        fo = open(fno, 'w')
        fo.write(new_string)

    def sweepOption(self, name, prefix, ext):
        fi = open(self.fn, 'r')
        string = fi.read()
        fi.close()
        while (self.opts[name].advance()):
            fno = '{}_{}.{}.{}'.format(prefix, name, self.opts[name].value, ext)
            self.writeCurrentState(fno)

    def sweepAllOptions(self, prefix, ext):
        for name in self.opts.keys():
            self.sweepOption(name, prefix, ext)
            
if __name__ == '__main__':
    opts = [
        Option('NLJ'     , '_NLJ_'     , [0, 1]        , 1) ,
        Option('LJA'     , '_LJA_'     , [0, 1]        , 1) ,
        Option('LJB'     ,  '_LJB_'    , [0, 1]        , 1) ,
        Option('NTORS'   , '_NTORS_'   , [0, 1]        , 1) ,
        Option('TORSNT'  , '_TORSNT_'  , [0,1]         , 1) ,
        Option('WTEMP'   , '_WTEMP_'   , [0, 300]      , 0) ,
        Option('DIST'    , '_DIST_'    , [1, 2, 3, 4]  , 4) ,
        Option('SELTYPE' , '_SEL_'     , [1, 2, 3, 4]  , 4) ,
        Option('CRTYPE'  , '_CR_'      , [1, 2]        , 1) ,
        Option('TORSTYPE', '_TORSTYPE_', [1, 2]        , 1) ,
        Option('MALG'    , '_MALG_'    , [0,1,2]       , 0)]
    
    inp = FileSweeper('input-template.inp', opts)
    inp.sweepAllOptions('output', 'inp')


