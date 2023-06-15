import unittest

from dsr_shelx import dbfile
from dsr_shelx.atomhandling import FindAtoms, NumberScheme
from dsr_shelx.resfile import ResList, ResListEdit


class TestReslist(unittest.TestCase):
    def setUp(self):
        # print ("setUp reslistTest executed!")
        self.atom1 = 'F10   4   -0.362398   0.278516   0.447770  11.00000   0.02302   0.04023 =\n'
        self.atom2 = '0.02897   0.00131  -0.01216   0.00374\n'
        self.res_file = 'tests/p21c.res'
        self.res_list = ResList(self.res_file)
        self.reslist = self.res_list.get_res_list()
        self.fa = FindAtoms(self.reslist)
        self.resi_str = 'RESI TOL 3'
        self.resi_list = ['RESI', 'TOL', '3']
        self.resi_list_class = ['RESI', 'TOL']
        self.resi_list_num = ['RESI', '3']

    def testrun_pos(self):
        self.assertEqual(self.fa.is_atom(self.atom1), ['F10', '4', '-0.362398', '0.278516', '0.447770'])
        self.assertEqual(self.fa.is_atom(self.atom2), [])

    def testrun_resinum(self):
        # print(self.fa.get_resi_definition_dict(self.resi_str))
        self.assertEqual(self.fa.get_resi_definition_dict(self.resi_str)['number'], '3')
        self.assertEqual(self.fa.get_resi_definition_dict(self.resi_list)['number'], '3')
        self.assertIsNone(self.fa.get_resi_definition_dict(self.resi_list_class)['number'])
        self.assertIs(self.fa.get_resi_definition_dict(self.resi_list_num)['number'], '3')

    def testrun_resinum_class(self):
        self.assertEqual(self.fa.get_resi_definition_dict(self.resi_str)['class'], 'TOL')
        self.assertEqual(self.fa.get_resi_definition_dict(self.resi_list)['class'], 'TOL')
        self.assertIs(self.fa.get_resi_definition_dict(self.resi_list_class)['class'], 'TOL')
        self.assertIsNone(self.fa.get_resi_definition_dict(self.resi_list_num)['class'])

    def testrun_get_atoms_resiclass(self):
        self.assertEqual(self.fa.get_atoms_resiclass('C1'), None)
        self.assertEqual(self.fa.get_atoms_resiclass('C1_1'), 'CF3')
        self.assertNotEqual(self.fa.get_atoms_resiclass('C1_1'), '1')

    def testrun_get_atoms_resinumber(self):
        self.assertEqual(self.fa.get_atoms_resinumber('C1_1'), '1')
        self.assertEqual(self.fa.get_atoms_resinumber('C1_1b'), '1')
        self.assertEqual(self.fa.get_atoms_resinumber('C1_b'), '0')
        self.assertEqual(self.fa.get_atoms_resinumber('C1_'), '0')
        self.assertNotEqual(self.fa.get_atoms_resinumber('C1_1'), 1)
        self.assertNotEqual(self.fa.get_atoms_resinumber('C1'), 0)
        self.assertEqual(self.fa.get_atoms_resinumber('C1_0'), '0')
        self.assertEqual(self.fa.get_atoms_resinumber('C1'), '0')


class remove_hydrogenTest(unittest.TestCase):
    def setUp(self):
        self.res_file = './collect_resi.res'
        self.res_list = ResList(self.res_file)
        self.reslist = self.res_list.get_res_list()
        self.fa = FindAtoms(self.reslist)

    def notworking_testrun_remove_adjacent_hydrogens(self):
        # self.res_list = ResList(self.res_file)
        # self.reslist =  self.res_list.get_res_list()
        # self.fa = FindAtoms(self.reslist)
        line1 = 'C1    1    0.462797    0.766414    0.415951    21.00000    0.01038    0.01899 =\n'
        line2 = '\n'
        sfac_table = ['C', 'O', 'F', 'AL', 'F', 'H']
        # print(self.reslist[37])
        # def showline(line):
        #    for num, i in enumerate(self.reslist):
        #        if num == line:
        #            print((i.strip('\n')))
        self.fa.remove_adjacent_hydrogens(['O1_1', 'C1_1'], sfac_table)
        # self.fa = FindAtoms(self.reslist)
        # showline(39)
        self.assertEqual(self.reslist[40], line1)
        self.assertEqual(self.reslist[41], line2)


class NumberSchemeTest(unittest.TestCase):
    def setUp(self):
        self.numbers = ['O1A', 'C1A', 'C2A', 'F1A', 'F2A', 'F3A', 'C3A',
                        'F4A', 'F5A', 'F6A', 'C4A', 'F7A', 'F8A', 'F9A']
        res_file = './tests/p21c.res'

        class Dsrp():
            resiflag = False

        dsrp = Dsrp()
        rl = ResList(res_file)
        reslist = rl.get_res_list()
        gdb = dbfile.ParseDB('./src/dsr_shelx/dsr_db.txt')
        fragment = 'OC(cf3)3'
        dbatoms = gdb.get_atoms(fragment)
        self.num = NumberScheme(reslist, dbatoms, dsrp)

    def testrun_get_numberscheme(self):
        numberscheme = self.num.get_fragment_number_scheme()
        self.assertListEqual(numberscheme, self.numbers)


class ResListEditTest(unittest.TestCase):
    def setUp(self):
        self.res_file = './tests/p21c.res'
        self.rl = ResList(self.res_file)
        self.res_list = self.rl.get_res_list()
        self.fa = FindAtoms(self.res_list)
        self.rle = ResListEdit(self.res_list, self.fa)

    def testrun_find_fvarlines(self):
        lines = self.rle.find_fvarlines()
        self.assertListEqual([17], lines)
        self.assertNotEqual(None, lines)
        self.assertNotEqual(False, lines)
