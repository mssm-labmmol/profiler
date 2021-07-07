from mock import patch
from numpy import full
import glob
import os
import unittest
import numpy as np
from ..profilerTools.linear_regression import LinearParameterOptimizer, LLS_SC
from ..profilerTools.profilerOpt.profilerOpt import ProfilerOptRunner

# =========== Patch all random functions (there are a lot). ==========

def patched_standard_normal(size=None):
    if (size is None):
        return 0.5
    elif (type(size) is int):
        return full(shape=(size,), fill_value=0.5)
    else:
        return full(shape=size, fill_value=0.5)
    
def patched_randn(*args):
    if len(args) == 0:
        return patched_standard_normal()
    else:
        return patched_standard_normal(size=args)

def patched_rand(*args):
    if len(args) == 0:
        return 0.5
    else:
        return full(shape=args, fill_value=0.5)

def patched_uniform(low, high, size=None):
    if (size is None):
        return 0.5 * (low + high)
    else:
        return full(shape=size, fill_value=0.5 * (low + high))
    
# numpy.standard_normal
patch_np_standard_normal = patch('numpy.random.standard_normal',
                                 side_effect=patched_standard_normal,
                                 spec=True)
patch_np_standard_normal.start()

# numpy.randn
patch_np_randn = patch('numpy.random.randn', side_effect=patched_randn,
                       spec=True)
patch_np_randn.start()

# numpy.randint
patch_np_randint = patch('numpy.random.randint',
                         side_effect=lambda low, high: low,
                         spec=True)
patch_np_randint.start()

# numpy.rand
patch_np_rand = patch('numpy.random.rand', side_effect=patched_rand,
                      spec=True)
patch_np_rand.start()

# numpy.uniform
patch_np_uniform = patch('numpy.random.uniform',
                         side_effect=patched_uniform, spec=True)
patch_np_uniform.start()

# random.uniform
patch_random_uniform = patch('random.uniform',
                             side_effect=patched_uniform, spec=True)
patch_random_uniform.start()

# random.gauss
patch_random_gauss = patch('random.gauss', side_effect=lambda mean,
                           stdev: mean + stdev *
                           patched_standard_normal(), spec=True)
patch_random_gauss.start()

# random.lognormvariate
patch_random_lognormvariate = patch('random.lognormvariate',
                                    side_effect=lambda mean, stdev:
                                    mean + stdev *
                                    patched_standard_normal(),
                                    spec=True)
patch_random_lognormvariate.start()

# random.choice
patch_random_choice = patch('random.choice', side_effect=lambda x:
                            x[0], spec=True)
patch_random_choice.start()

# random.random
patch_random_random = patch('random.random', return_value=0.5)
patch_random_random.start()

# random.sample
patch_random_sample = patch('random.sample', side_effect=lambda x, k:
                            x[:k], spec=True)

patch_random_sample.start()

# random.randint
patch_random_randint = patch('random.randint', side_effect=lambda low,
                             high: low, spec=True)
patch_random_randint.start()

# random.randrange
patch_random_randrange = patch('random.randrange', side_effect=lambda
                               low, high: low, spec=True)
patch_random_randrange.start()

# random.shuffle
patch_random_shuffle = patch('random.shuffle', side_effect=lambda x:
                             x, spec=True)

patch_random_shuffle.start()


class MockMultiprofile:


    def __init__(self):
        pass


    def calculateNumberOfConfigurations(self):
        return [2,1,3]

    def getNonoptEnergy(self):
        return np.array([0,0,0])


class TestLPO(unittest.TestCase):

    def test_build_S(self):
        Sref = np.array([
            [1,0,0],
            [1,0,0],
            [0,1,0],
            [0,0,1],
            [0,0,1],
            [0,0,1],
        ])
        lpo = LinearParameterOptimizer(MockMultiprofile())
        np.testing.assert_array_equal(lpo.Sm, Sref)


    def test_fit(self):
        args = [
            '-c', 'tests/data/profilerOpt_bu.xyz', 'tests/data/profilerOpt_pe.xyz',
            '-r', 'tests/data/profilerOpt_bu.dat', 'tests/data/profilerOpt_pe.dat',
            '-t', 'tests/data/profilerOpt_bu.stp', 'tests/data/profilerOpt_pe.stp',
            '-i', 'tests/data/profilerOpt_run_short.inp',
            '-op', 'tests/output_test_linear_regression',
        ]


        # run job
        job = ProfilerOptRunner(args)


        # ======================= Test matrix A ======================

        lpo = LinearParameterOptimizer(job.get_optimal_individual())


        Am = lpo._build_matrix_A()

        # slice to remove energy-shift columns
        Aref = np.loadtxt('tests/data/linear_regression_X.dat')[:,:-2]

        # compare
        #np.testing.assert_array_almost_equal(Am, Aref, decimal=3)

        # ========================= Test fit =========================
        
        # test fitting
        lpo.fit(job.get_target_data())
        pars = lpo.model.coef_[0,:-2]
        pars_ref = [6.558520e-04, 1.973677e-06, 1.254371e-03, 1.981070e-06, 6.955627e+00]

        #np.testing.assert_array_almost_equal(pars, pars_ref, decimal=3)

        # =================== Test self-consistent ===================

        lls_sc = LLS_SC(job.get_optimal_individual())
        lls_sc.run(40, 1.0e-04, job.get_target_data())
        
        # remove output files
        for fn in glob.glob('tests/output_test_linear_regression*'):
            os.remove(fn)
        
        
if __name__ == '__main__':
    unittest.main()
    patch_np_standard_normal.stop()
    patch_np_randn.stop()
    patch_np_randint.stop()
    patch_np_rand.stop()
    patch_np_uniform.stop()
    patch_random_uniform.stop()
    patch_random_gauss.stop()
    patch_random_lognormvariate.stop()
    patch_random_choice.stop()
    patch_random_random.stop()
    patch_random_sample.stop()
    patch_random_randint.stop()
    patch_random_randrange.stop()
    patch_random_shuffle.stop()
