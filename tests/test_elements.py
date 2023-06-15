import unittest

from dsr_shelx.atomhandling import get_atomtypes, SfacTable, Elem_2_Sfac
from dsr_shelx.atoms import Element
from dsr_shelx.resfile import ResList


class TestElement(unittest.TestCase):
    def setUp(self):
        self.element = Element()

    def test_get_atomic_number_valid(self):
        self.assertEqual(self.element.get_atomic_number('F'), 9)
        self.assertEqual(self.element.get_atomic_number('C'), 6)

    def test_get_number_from_q_peak(self):
        self.assertEqual(self.element.get_atomic_number('Q'), 6)

    def test_get_number_from_charged_atom(self):
        self.assertEqual(None, self.element.get_atomic_number('O1-'))

if __name__ == '__main__':
    unittest.main()


class ElementsTest(unittest.TestCase):
    def setUp(self):
        self.el = Element()

    def testrun_get_atomic_number(self):
        num = self.el.get_atomic_number('Fe')
        ele = self.el.get_element(26)
        self.assertEqual(num, 26)
        self.assertEqual(ele, 'Fe')


class SfacTableTest(unittest.TestCase):
    def setUp(self):
        self.res_file = './tests/collect_resi.res'
        self.res_list = ResList(self.res_file)
        self.reslist = self.res_list.get_res_list()
        self.fragment_atom_types = ['C', 'H', 'N', 'Na']
        self.reference = ['C', 'O', 'F', 'AL', 'FE', 'H', 'N', 'NA']
        self.bad_atomtypes = ['C', 'H', 'Naa', 'XY']

    def testrun_SFacTable(self):
        sf = SfacTable(self.reslist, self.fragment_atom_types)
        sfac_list = sf.set_sfac_table()
        self.assertListEqual(self.reference, sfac_list)
        sfbad = SfacTable(self.reslist, self.bad_atomtypes)
        with self.assertRaises(SystemExit):
            sfbad.set_sfac_table()


class Elem_2_SfacTest(unittest.TestCase):
    def setUp(self):
        self.sfac_table = ['C', 'H', 'N', 'Na']
        self.testlist = ['C', 'H', 'N', 'Na', 'NA', 'n']
        self.upper_testlist = ['C', 'H', 'N', 'NA', 'NA', 'N']
        self.returnlist = [1, 2, 3, 4, 4, 3]

    def testrun_elem_2_sfac(self):
        e2s = Elem_2_Sfac(self.sfac_table)
        outcome = []
        for i in self.testlist:
            num = e2s.elem_2_sfac(i)
            outcome.append(num)
        self.assertListEqual(outcome, self.returnlist)
        self.assertEqual(e2s.elem_2_sfac('Caa'), None)

    def testrun_sfac_2_elem(self):
        e2s = Elem_2_Sfac(self.sfac_table)
        outcome = []
        for i in self.returnlist:
            element = e2s.sfac_2_elem(i)
            outcome.append(element)
        self.assertListEqual(outcome, self.upper_testlist)


class get_atomtypesTest(unittest.TestCase):
    def setUp(self):
        self.dbatoms = [['O1', 3, '-0.01453', '1.66590', '0.10966'],
                        ['C1', 1, '-0.00146', '0.26814', '0.06351']]
        self.bad_dbatoms1 = [['lO1', 3, '-0.01453', '1.66590', '0.10966'],
                             ['C1', 1, '-0.00146', '0.26814', '0.06351']]
        self.bad_dbatoms2 = [['O1', 3, '-0.01453', '1.66590', '0.10966'],
                             ['1', 1, '-0.00146', '0.26814', '0.06351']]

    def testrun_pos(self):
        self.assertEqual(get_atomtypes(self.dbatoms), ['O', 'C'])

    def testrun_neg1(self):
        with self.assertRaises(KeyError):
            get_atomtypes(self.bad_dbatoms1)

    def testrun_neg2(self):
        with self.assertRaises(KeyError):
            get_atomtypes(self.bad_dbatoms2)
