import unittest

from dsr_shelx.atoms import atoms
from dsr_shelx.dbfile import invert_atomic_coordinates


class invert_fragmentTest(unittest.TestCase):
    def setUp(self):
        self.dbatoms = [['O1', 3, -0.01453, 1.66590, 0.10966], ['C1', 1, -0.00146, 0.26814, 0.06351],
                        ['C2', 1, -1.13341, -0.23247, -0.90730], ['F1', 4, -2.34661, -0.11273, -0.34544]]
        self.dbatoms2 = [['O1', 3, -0.01453, 1.66590, 0.10966], ['C1', 1, -0.00146, 0.26814, 0.06351],
                         ['C2', 1, -1.13341, -0.23247, -0.90730], ['F1', 4, -2.34661, -0.11273, -0.34544]]
        self.inv_dbatoms = [['O1', 3, 0.01453, -1.6659, -0.10966], ['C1', 1, 0.00146, -0.26814, -0.06351],
                            ['C2', 1, 1.13341, 0.23247, 0.9073], ['F1', 4, 2.34661, 0.11273, 0.34544]]

    def testrun_invert_dbatoms_coordinates(self):
        inverted = invert_atomic_coordinates(self.dbatoms)
        self.assertListEqual(self.inv_dbatoms, inverted)
        self.assertNotEqual(self.dbatoms2[0][2], inverted[0][2])


class atomsTest(unittest.TestCase):

    def testrun_num_of_atoms(self):
        num_of_atoms = len(atoms)
        self.assertEqual(num_of_atoms, 93)
        self.assertNotEqual(num_of_atoms, 42)
