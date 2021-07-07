from mock import patch
from numpy import full
import glob
import os
import unittest
import numpy as np
from ..profilerTools.readopts import Global_Type_IndexConverter


class TestGlobal_Type_IndexConverter(unittest.TestCase):

    def test(self):

        nTors = 3
        nLJ   = 2
        kMask = [
            [1, 1, 0, 0, 0, 1],
            [0, 0, 1, 0, 0, 0],
            [0, 1, 0, 0, 1, 0],
        ]
        phiMask = [
            [1, 1, 0, 0, 0, 1],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0],
        ]
        LJMask = [
            [1, 0],
            [1, 1],
        ]


        gb = Global_Type_IndexConverter(nTors, nLJ, kMask, phiMask, LJMask)

        self.assertEqual(gb.global_to_type(0), (0, 'cs6'))
        self.assertEqual(gb.global_to_type(1), (1, 'cs6'))
        self.assertEqual(gb.global_to_type(6), (2, 'phi_1'))
        self.assertEqual(gb.global_to_type(7), (2, 'phi_2'))
        self.assertEqual(gb.global_to_type(8), (2, 'phi_6'))
        self.assertEqual(gb.global_to_type(9), (3, 'k_3'))
        self.assertEqual(gb.type_to_global(1), [1,2])
        self.assertEqual(gb.type_to_global(2), [3,4,5,6,7,8])
        self.assertEqual(gb.type_to_global(3), [9])
        self.assertEqual(gb.get_k_phi_pairs(), [(3,6), (4,7), (5,8), (11,12)])
        
        
if __name__ == '__main__':
    unittest.main()
