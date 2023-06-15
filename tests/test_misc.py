import unittest

from dsr_shelx import atomhandling, misc
from dsr_shelx.restraints import format_atom_names
from tests.test_all import db_testhead, wraphead, atomnames, formated_atoms, part_atoms, resi_atoms, resi_atoms2, cell, \
    coord1, coord2


class MiscTest(unittest.TestCase):
    def setUp(self):
        self.dbatoms = [['O1', '3', '-0.01453', '1.66590', '0.10966'],
                        ['C1', '1', '-0.00146', '0.26814', '0.06351'],
                        ['C2', '1', '-1.13341', '-0.23247', '-0.90730'],
                        ['Sn1', '4', '-2.34661', '-0.11273', '-0.34544']]
        self.strdbatoms = ['O1 3 -0.01453 1.66590 0.10966 11.00 0.05',
                           'C1 1 -0.00146 0.26814 0.06351 11.00 0.05',
                           'C2 1 -1.13341 -0.23247 -0.90730 11.00 0.05',
                           'Sn1 4 -2.34661 -0.11273 -0.34544 -21.00 0.05']
        self.residue_atoms = ['O1 3 -0.01453 1.66590 0.10966 11.00 0.05',
                              'C1 1 -0.00146 0.26814 0.06351 11.00 0.05',
                              'RESI 4 BENZ',
                              'C2 1 -1.13341 -0.23247 -0.90730 11.00 0.05',
                              'Sn1 4 -2.34661 -0.11273 -0.34544 -21.00 0.05']
        self.string = 'O1   3   -0.01453   1.66590   0.10966\nC1   1   -0.00146   0.26814   0.06351\nC2   1   -1.13341   -0.23247   -0.90730\nSn1   4   -2.34661   -0.11273   -0.34544'
        self.multi = 'rem dsr put oc(cf3)3 with o1 c1 c2 c3 c4 on O1_3 c1_3 q6 Q4 q7 resi cf3 =\
                        PART 2 occ -31'

    def testrun_find_line_of_residue(self):
        self.maxDiff = None
        number, string = line = misc.find_line_of_residue(self.residue_atoms, '4')
        self.assertNotEqual(number, 3)
        self.assertEqual(number, 2)
        self.assertEqual(len(line), 2)
        self.assertEqual(string, 'RESI 4 BENZ')

    def testrun_get_atoms(self):
        self.maxDiff = None
        noatoms = atomhandling.get_atoms([])
        self.assertEqual(noatoms, [])
        # atoms = atomhandling.get_atoms(self.dbatoms)
        # self.assertListEqual(atoms, self.dbatoms)
        stratoms = atomhandling.get_atoms(self.strdbatoms)
        self.assertListEqual(stratoms, self.dbatoms)
        tst = atomhandling.get_atoms(['O1 3 -0.01453 1.66590 0.10966'])
        tst2 = atomhandling.get_atoms(['O1 3 -0.01453 1.66590'])
        tst3 = atomhandling.get_atoms('O1 3 -0.01453 1.66590  0.10966')
        self.assertEqual(tst, [['O1', '3', '-0.01453', '1.66590', '0.10966']])
        self.assertEqual(tst2, [])
        self.assertEqual(tst3, [])

    def testrun_ll_to_string(self):
        llstr = misc.ll_to_string(self.dbatoms)
        self.assertEqual(llstr, self.string)

    def testrun_multiline(self):
        self.maxDiff = None
        mult = misc.multiline_test(self.multi)
        mult2 = misc.multiline_test('O1 3 -0.01453 1.66590 0.10966 11.00 0.05')
        with self.assertRaises(AttributeError):
            misc.multiline_test([['O1', '3', '-0.01453', '1.66590', '0.10966']])
        self.assertEqual(mult, True)
        self.assertEqual(mult2, False)

    def testrun_find_multi_lines(self):
        self.maxDiff = None
        found = misc.find_multi_lines(self.strdbatoms, '[A-z][0-9]')
        self.assertEqual(found, [0, 1, 2])
        found2 = misc.find_multi_lines(self.strdbatoms, '[A-z][0-9] khgf')
        self.assertEqual(found2, [])
        with self.assertRaises(TypeError):
            misc.find_multi_lines(self.dbatoms, '[A-z][0-9]')

    def testrun_wrap_headlines(self):
        self.maxDiff = None
        unwraped = ['SADI C1 C2 C1 C3 C1 C4',
                    'SADI F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4',
                    'SADI 0.04 C2 C3 C3 C4 C2 C4', 'SADI 0.04 O1 C2 O1 C3 O1 C4',
                    'SADI 0.04 F1 F2 F2 F3 F3 F1 F4 F5 F5 F6 F6 F4 F7 F8 F8 F9 F9 F7',
                    'SADI 0.04 F1 C1 F2 C1 F3 C1 F4 C1 F5 C1 F6 C1 F7 C1 F8 C1 F9 C1',
                    'SIMU O1 > F9', 'RIGU O1 > F9']
        head = misc.wrap_headlines(db_testhead)
        self.assertListEqual(head, wraphead)
        unwrap = misc.unwrap_head_lines(wraphead)
        self.assertListEqual(unwraped, unwrap)

    def disabled_testrun_which(self):
        which = misc.which('notepad')[0]
        self.assertEqual(which, 'C:\\Windows\\system32\\notepad.exe')

    def testrun_zero(self):
        zero = misc.zero(3, 3)
        zer0 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        self.assertListEqual(zero, zer0)

    def testrun_format_atom_names(self):
        names = format_atom_names(atomnames, part=2, resinum=4)
        names2 = format_atom_names(atomnames, part='2', resinum='4')
        names3 = format_atom_names(atomnames, part='2')
        names4 = format_atom_names(atomnames, resinum=91)
        names5 = format_atom_names(atomnames, part='2a')
        self.assertListEqual(formated_atoms, names)
        self.assertListEqual(formated_atoms, names2)
        self.assertListEqual(part_atoms, names3)
        self.assertListEqual(resi_atoms, names4)
        self.assertListEqual(atomnames, names5)

    def testrun_remove_partsymbol(self):
        def removep(at):
            atoms = []
            for i in at:
                atoms.append(misc.remove_partsymbol(i))
            return atoms

        remove1 = removep(formated_atoms)
        remove2 = removep(part_atoms)
        remove3 = removep(resi_atoms)
        self.assertListEqual(resi_atoms2, remove1)
        self.assertListEqual(atomnames, remove2)
        self.assertListEqual(resi_atoms, remove3)

    def testrun_atomic_distance(self):
        cellf = [float(i) for i in cell]
        dst = misc.atomic_distance(coord1, coord2, cellf)
        self.assertAlmostEqual(6.05249787959, dst, 5)

    def testrun_frac_to_cart(self):
        coord1 = (-0.186843, 0.282708, 0.526803)
        cellf = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
        cart = (-2.741505423999065, 5.909586678000002, 10.775200700893734)
        N1 = misc.frac_to_cart(coord1, cellf)
        for x, y in zip(cart, N1):
            self.assertAlmostEqual(x, y, 8)

    def testrun_subtract(self):
        m1 = [2, 0, 0]
        m2 = [1, 1.5, 0.104]
        sub = misc.subtract_vect(m1, m2)
        self.assertEqual(sub, (1, -1.5, -0.104))

    def testrun_vol_tetrahedron(self):
        # cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
        a = (0.838817, 0.484526, 0.190081)  # a ist um 0.01 ausgelenkt
        b = (0.875251, 0.478410, 0.256955)
        c = (0.789290, 0.456520, 0.301616)
        d = (0.674054, 0.430194, 0.280727)
        vol = misc.vol_tetrahedron(a, b, c, d, cell)
        self.assertAlmostEqual(0.063352, vol, 5)

    def testrun_dice_coefficient(self):
        string = 'Die Katze lauft im Schnee! CF3, Benzene'
        srch1 = 'Katz'
        srch2 = 'Benze'
        srch3 = 'Hallo du'
        srch4 = '124'
        dice1 = misc.dice_coefficient(string, srch1)
        dice2 = misc.dice_coefficient(string, srch2)
        dice3 = misc.dice_coefficient(string, srch3)
        dice4 = misc.dice_coefficient(string, srch4)
        self.assertEqual(dice1, 0.162162)
        self.assertEqual(dice2, 0.210526)
        self.assertEqual(dice3, 0.0)
        self.assertEqual(dice4, 0.0)
