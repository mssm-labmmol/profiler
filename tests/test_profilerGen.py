from mock import patch
from numpy import full
import os
import glob
import unittest
from   numpy import loadtxt
from ..profilerGen.profilerGen import ProfilerGenRunner

class TestProfilerGen(unittest.TestCase):

    def test_integrated(self):
        args = [
            '-c', 'tests/data/profilerGen_pe.xyz',
            '-t', 'tests/data/profilerGen_pe.stp',
            '-dk', '5000', '5000',
            '-op', 'tests/output_test_profilerGen',
            '-dr', '-180', '90', '180',
            '-dr', '-180', '120', '180',
            '-min', '1',
            '-dx0', '0.05',
            '-dxm', '0.20',
            '-dele', '1.0e-06',
            '-nsteps', '10000',
        ]
        refdata = list(loadtxt('tests/data/output_test_profilerGen.dat',
                          usecols=(0,)))
        

        # run job
        job = ProfilerGenRunner(args)

        # remove output files
        for fn in glob.glob('tests/output_test_profilerGen*'):
            os.remove(fn)

        # compare
        for calc, ref in zip(job.get_energies(), refdata):
            self.assertAlmostEqual(calc, ref, places=4)
        
if __name__ == '__main__':
    unittest.main()
