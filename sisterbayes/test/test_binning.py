#! /usr/bin/env python

import unittest
import bisect
from sisterbayes import utility

class BisectTestCase(unittest.TestCase):

    def validateBinIndex(self, value, bin_size):
        obs = utility.bin_index(value, bin_size)
        factor = 1e6
        check_array = [i/float(factor) for i in range(0, int(value * 2 * factor), int(bin_size * factor))]
        self.assertEqual(bisect.bisect_left(check_array, value), obs)

    def test_int(self):
        for value in [3, 4, 7, 9, 17, 1, 0, 15, 2012, 23, 1, 83]:
            for bin_size in [1,2,3,4,5,7,11,13,100]:
                self.validateBinIndex(value, bin_size)

    def test_float(self):
        for value in [3.14, 0.04, 0.12, 9.112, 0.0017, 0.00511, 0.12, 15.173741, 2.18182, 0.123, 0.101, 0.00283]:
            for bin_size in [0.1, 0.001, 0.5, 1.0, 0.015, 0.13, 0.00001]:
                self.validateBinIndex(value, bin_size)

if __name__ == "__main__":
    unittest.main()
