#!/usr/bin/env python3

# The test relies on the installed version of gslpyridge.

from gslpyridge import fit
import unittest
import numpy as np


class TestGSLRidge(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.X_real = np.array([
            [1.0, 0.4, 1.7],
            [0.2, 2.2, 0.8],
            [0.1, 1.2, 3.2],
            [4.2, 1.1, 0.2],
        ])
        cls.p_real = np.array([0.30, 0.15, 0.25]).reshape(-1, 1)
        cls.Y_real = np.matmul(cls.X_real, cls.p_real)

        
    def test_no_reg(self):
        no_reg_y = fit(X=self.X_real, y=self.Y_real, lamb=0)
        np.testing.assert_almost_equal(self.p_real, no_reg_y, decimal=3)


if __name__ == '__main__':
    unittest.main()
