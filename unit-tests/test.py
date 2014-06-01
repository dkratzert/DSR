#/usr/bin/env python
#-*- encoding: utf-8 -*-


import sys, os
import unittest
from ..atomhandling import get_atomtypes

dbatoms = [['O1', 3, '-0.01453', '1.66590', '0.10966'], ['C1', 1, '-0.00146', '0.26814', '0.06351']]

class get_atomtypesTest(unittest.TestCase):

    def testrun(self):
        self.assertEqual(get_atomtypes(dbatoms), ['O', 'C'])

#if __name__ == "__main__":
#    unittest.main()

if __name__ == "__main__" and __package__ is None:
    __package__ = "test.py"
    