from mock import patch
from numpy import full
import glob
import os
import unittest
import numpy as np
from ..profilerTools.geometry_factor import *
from ..profilerTools.configuration import ensemble


class Test_geometry_factor(unittest.TestCase):

    def test(self):
        ens = ensemble()

        ens.readFromTrajectory('tests/data/profilerGen_pe.xyz')
        conf = ens[0]

        selstring_1 = "cs12"
        atom_idxs_1 = np.array([[2,1], [0,4]], dtype='int32')

        selstring_2 = "phi_2"
        atom_idxs_2 = np.array([[2,0,3,1], [4,1,3,0]], dtype='int32')

        dofs_1 = calculate_dofs(conf, atom_idxs_1, selstring_1)
        dofs_2 = calculate_dofs(conf, atom_idxs_2, selstring_2)

        # check
        np.testing.assert_array_almost_equal(dofs_1, [0.2939545, 0.39464652], decimal=3)
        np.testing.assert_array_almost_equal(dofs_2, [0.0, 180.0], decimal=3)

        
if __name__ == '__main__':
    unittest.main()
