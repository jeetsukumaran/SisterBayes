#! /usr/bin/env python

import unittest
import bisect
from sisterbayes import utility

class BisectTestCase(unittest.TestCase):

    def validateBinIndex(self, value, bin_size):
        obs = utility.bin_index(value, bin_size)
        check_array = [i for i in range(0, value * 2, bin_size)]
        self.assertEqual(bisect.bisect_left(check_array, value), obs)

    def test1(self):
        for value in [3, 4, 7, 9, 17, 1, 0, 15, 2012, 23, 1, 83]:
            for bin_size in [1,2,3,4,5,7,11,13,100]:
                self.validateBinIndex(value, bin_size)

if __name__ == "__main__":
    unittest.main()
