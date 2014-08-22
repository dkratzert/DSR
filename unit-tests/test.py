#/usr/bin/env python
#-*- encoding: utf-8 -*-


# test Afix: write_dbhead_to_file()
# test resi module

import unittest
from dsr import VERSION
from afix import InsertAfix
from atomhandling import get_atomtypes, FindAtoms, check_source_target, \
    rename_dbhead_atoms, SfacTable, Elem_2_Sfac, NumberScheme
from atoms import Element
from dbfile import global_DB, invert_dbatoms_coordinates, ReadDB, ImportGRADE
from dsrparse import DSR_Parser
import misc
from resfile import ResList, ResListEdit
from resi import Resi
from export import Export

db_testhead = ['SADI C1 C2 C1 C3 C1 C4',
               'SADI F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4',
               'SADI 0.04 C2 C3 C3 C4 C2 C4', 'SADI 0.04 O1 C2 O1 C3 O1 C4',
               'SADI 0.04 F1 F2 F2 F3 F3 F1 F4 F5 F5 F6 F6 F4 F7 F8 F8 F9 F9 F7',
               'SADI 0.1 F1 C1 F2 C1 F3 C1 F4 C1 F5 C1 F6 C1 F7 C1 F8 C1 F9 C1',
               'SIMU O1 > F9', 'RIGU O1 > F9']

wraphead = ['SADI C1 C2 C1 C3 C1 C4\n',
         'SADI F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4 F4 C3 F5 C3 F6 C3 =\n   F7 C4 F8 C4 F9 C4\n',
         'SADI 0.04 C2 C3 C3 C4 C2 C4\n', 'SADI 0.04 O1 C2 O1 C3 O1 C4\n',
         'SADI 0.04 F1 F2 F2 F3 F3 F1 F4 F5 F5 F6 F6 F4 F7 F8 F8 F9 F9 F7\n',
         'SADI 0.1 F1 C1 F2 C1 F3 C1 F4 C1 F5 C1 F6 C1 F7 C1 F8 C1 F9 C1\n',
         'SIMU O1 > F9\n', 'RIGU O1 > F9\n']

atomnames = ['O1', 'C1', 'C2', 'F1', 'F2', 'F3', 'C3', 'F4', 'F5', 'F6', 'C4', 'F7', 'F8', 'F9']
formated_atoms = ['O1_4B', 'C1_4B', 'C2_4B', 'F1_4B', 'F2_4B', 'F3_4B', 'C3_4B', 'F4_4B', 'F5_4B', 'F6_4B', 'C4_4B', 'F7_4B', 'F8_4B', 'F9_4B']
part_atoms = ['O1_B', 'C1_B', 'C2_B', 'F1_B', 'F2_B', 'F3_B', 'C3_B', 'F4_B', 'F5_B', 'F6_B', 'C4_B', 'F7_B', 'F8_B', 'F9_B']
resi_atoms = ['O1_91', 'C1_91', 'C2_91', 'F1_91', 'F2_91', 'F3_91', 'C3_91', 'F4_91', 'F5_91', 'F6_91', 'C4_91', 'F7_91', 'F8_91', 'F9_91']
resi_atoms2 = ['O1_4', 'C1_4', 'C2_4', 'F1_4', 'F2_4', 'F3_4', 'C3_4', 'F4_4', 'F5_4', 'F6_4', 'C4_4', 'F7_4', 'F8_4', 'F9_4']

coord1 = (0.100759, 0.449528, 0.438703) #F4
coord1a = ('0.100759', '0.449528', '0.438703') #F4
coord2 = (-0.362398, 0.278516, 0.447770) #F10 6.052A
cell = [10.5086, 20.9035, 20.5072, 90, 94.13, 90]
cells = ['10.5086', '20.9035', '20.5072', '90', '94.13', '90']


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
        with self.assertRaises(SystemExit):
            get_atomtypes(self.bad_dbatoms1)
    def testrun_neg2(self):
        with self.assertRaises(SystemExit):
            get_atomtypes(self.bad_dbatoms2)

class reslistTest(unittest.TestCase):
    def setUp(self):
        #print ("setUp reslistTest executed!")
        self.atom1 = 'F10   4   -0.362398   0.278516   0.447770  11.00000   0.02302   0.04023 =\n'
        self.atom2 = '0.02897   0.00131  -0.01216   0.00374\n'
        self.res_file = 'p21c.res'
        self.res_list = ResList(self.res_file)
        self.reslist =  self.res_list.get_res_list()
        self.fa = FindAtoms(self.reslist)
        self.resi_str = 'RESI TOL 3'
        self.resi_list = ['RESI', 'TOL', '3']
        self.resi_list_class = ['RESI', 'TOL']
        self.resi_list_num = ['RESI', '3']

    def testrun_pos(self):
        self.assertEqual(self.fa.get_atom(self.atom1), ['F10', '4', '-0.362398', '0.278516', '0.447770'])
        self.assertIsNone(self.fa.get_atom(self.atom2))

    def testrun_resinum(self):
        #print(self.fa.get_resinum(self.resi_str))
        self.assertEqual(self.fa.get_resinum(self.resi_str)['number'], '3')
        self.assertEqual(self.fa.get_resinum(self.resi_list)['number'], '3')
        self.assertIsNone(self.fa.get_resinum(self.resi_list_class)['number'])
        self.assertIs(self.fa.get_resinum(self.resi_list_num)['number'], '3')

    def testrun_resinum_class(self):
        self.assertEqual(self.fa.get_resinum(self.resi_str)['class'], 'TOL')
        self.assertEqual(self.fa.get_resinum(self.resi_list)['class'], 'TOL')
        self.assertIs(self.fa.get_resinum(self.resi_list_class)['class'], 'TOL')
        self.assertIsNone(self.fa.get_resinum(self.resi_list_num)['class'])

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


class collect_residuesTest(unittest.TestCase):
    def setUp(self):
        self.res_file = './collect_resi.res'
        self.res_list = ResList(self.res_file)
        self.reslist =  self.res_list.get_res_list()
        self.fa = FindAtoms(self.reslist)

    def testrun_collect_residues(self):
        #print(self.fa.collect_residues())
        collected_resi = {'0': [['FE1', ['0.100001', '0.200002', '0.300003'], 19, None]],
                          '1': [['O1', ['0.584527', '0.749093', '0.406892'], 34, 'CCF3'],
                                ['H2', ['0.2', '0.3', '0.4'], 36, 'CCF3'],
                                ['C1', ['0.462797', '0.766414', '0.415951'], 37, 'CCF3'],
                                ['H1', ['0.2', '0.3', '0.4'], 40, 'CCF3']]}
        resi_null = [['FE1', ['0.100001', '0.200002', '0.300003'], 19, None]]
        resi_fe = 'FE1'
        self.assertEqual(self.fa.collect_residues(), collected_resi)
        self.assertEqual(self.fa.collect_residues()['0'], resi_null)
        self.assertEqual(self.fa.collect_residues()['0'][0][0], resi_fe)

    def testrun_get_atomcoordinates(self):
        atoms_fe = ['Fe1']
        atoms_fe_0 = ['Fe1_0']
        atoms = ['Fe1', 'C1_1', 'O1_1']
        fe_upper = 'FE1'
        fe = {'FE1': ['0.100001', '0.200002', '0.300003']}
        fe_float = {'FE1': [0.100001, 0.200002, 0.300003]}
        coords = ['0.100001', '0.200002', '0.300003']
        atoms_coords = {'C1_1': ['0.462797', '0.766414', '0.415951'],
                        'O1_1': ['0.584527', '0.749093', '0.406892'],
                        'FE1': ['0.100001', '0.200002', '0.300003']}

        self.assertEqual(self.fa.get_atomcoordinates(atoms_fe), fe)
        self.assertNotEqual(self.fa.get_atomcoordinates(atoms_fe), fe_float)
        self.assertEqual(list(self.fa.get_atomcoordinates(atoms_fe).keys())[0], fe_upper)
        self.assertEqual(self.fa.get_atomcoordinates(atoms_fe_0)['FE1_0'], coords)
        self.assertEqual(self.fa.get_atomcoordinates(atoms), atoms_coords)
        with self.assertRaises(KeyError):
            # raises KeyError, because FE1 is in residue 0
            self.fa.get_atomcoordinates(atoms_fe_0)['FE1_1']

    def testrun_get_atom_line_number(self):
        line_str = ['34']
        line_int = [34]
        self.assertNotEqual(self.fa.get_atom_line_numbers(['O1_1']), line_str)
        self.assertEqual(self.fa.get_atom_line_numbers(['O1_1']), line_int)
        self.assertEqual(self.fa.get_atom_line_numbers(['Fe1_0', 'C1_1']), [19, 37])
        self.assertEqual(self.fa.get_atom_line_numbers(['C1_1', 'Fe1_0']), [37, 19])

    def testrun_check_source_target(self):
        source = ['C1', 'C3', 'O1', 'F1']
        source_short = ['C1', 'C3', 'O1']
        source_fail = ['C15', 'C3', 'O1', 'F1']
        target = ['q1', 'C3', 'O1_2', 'F1']
        target_num = ['q1', 'C3', 'O1_2']
        dbatoms = [['O1', 3, '-0.01453', '1.66590', '0.10966'], ['C1', 1, '-0.00146', '0.26814', '0.06351'],
                   ['C2', 1, '-1.13341', '-0.23247', '-0.90730'], ['F1', 4, '-2.34661', '-0.11273', '-0.34544'],
                   ['F2', 4, '-0.96254', '-1.50665', '-1.29080'], ['F3', 4, '-1.12263', '0.55028', '-2.01763'],
                   ['C3', 1, '1.40566', '-0.23179', '-0.43131'], ['F4', 4, '2.38529', '0.42340', '0.20561'],
                   ['F5', 4, '1.53256', '0.03843', '-1.75538'], ['F6', 4, '1.57833', '-1.55153', '-0.25035'],
                   ['C4', 1, '-0.27813', '-0.21605', '1.52795'], ['F7', 4, '0.80602', '-0.03759', '2.30431'],
                   ['F8', 4, '-0.58910', '-1.52859', '1.53460'], ['F9', 4, '-1.29323', '0.46963', '2.06735']]
        self.assertEqual(check_source_target(source, target, dbatoms), True)
        with self.assertRaises(SystemExit):
            check_source_target(source, target_num, dbatoms)
        with self.assertRaises(SystemExit):
            check_source_target(source_fail, target, dbatoms)
        with self.assertRaises(SystemExit):
            check_source_target(source_short, target, dbatoms)


class remove_hydrogenTest(unittest.TestCase):
    #def setUp(self):
        #self.res_file = 'unit-tests/collect_resi.res'
        #self.res_list = ResList(self.res_file)
        #self.reslist =  self.res_list.get_res_list()
        #self.fa = FindAtoms(self.reslist)

    def notworking_testrun_remove_adjacent_hydrogens(self):
        self.res_file = 'unit-tests/collect_resi.res'
        self.res_list = ResList(self.res_file)
        self.reslist =  self.res_list.get_res_list()
        self.fa = FindAtoms(self.reslist)
        #line1 = 'C1    1    0.462797    0.766414    0.415951    21.00000    0.01038    0.01899 =\n'
        #line2 = '\n'
        sfac_table = ['C', 'O', 'F', 'AL', 'F', 'H']
        #print(self.reslist[37])
        def showline(line):
            for num, i in enumerate(self.reslist):
                if num == line:
                    print((i.strip('\n')))
        self.fa.remove_adjacent_hydrogens(['O1_1', 'C1_1'], sfac_table)
        #self.fa = FindAtoms(self.reslist)
        showline(39)
        #self.assertEqual(self.reslist[36], line1)
        #self.assertEqual(self.reslist[39], line2)


class rename_DBHeadatomsTest(unittest.TestCase):
    def setUp(self):
        self.head =['SADI_CF3 C1 C2 C1 C3 C1 C4', 'SADI_CF3 F1 C2 F2 C2 F3 C2 F4 C3'
               ' F5 C3 F6 C3 F7 C4 F8 C4 F9 C4', 'SADI_CF3 0.04 C2 C3 C3 C4 C2 C4',
               'SADI_CF3 0.04 O1 C2 O1 C3 O1 C4', 'SADI_CF3 0.04 F1 F2 F2 F3 F3'
               ' F1 F4 F5 F5 F6 F6 F4 F7 F8 F8 F9 F9 F7', 'SADI_CF3 0.1 F1 C1 F2 '
               'C1 F3 C1 F4 C1 F5 C1 F6 C1 F7 C1 F8 C1 F9 C1', 'SIMU_CF3 O1 > F9',
               'RIGU_CF3 O1 > F9', 'RESI 4 CF3']
        self.old = ('O1', 'C1', 'C2', 'F1', 'F2', 'F3', 'C3', 'F4', 'F5', 'F6', 'C4', 'F7', 'F8', 'F9')
        self.new = ('F9A', 'F8A', 'F7A', 'C4A', 'F6A', 'F5A', 'F4A', 'C3A', 'F3A', 'F2A', 'F1A', 'C2A', 'C1A', 'O1A')
        self.new_rev = tuple(reversed(['F9A', 'F8A', 'F7A', 'C4A', 'F6A', 'F5A', 'F4A', 'C3A', 'F3A', 'F2A', 'F1A', 'C2A', 'C1A', 'O1A']))

    def testrun_rename_dbheadatoms(self):
        newhead = rename_dbhead_atoms(self.new, self.old, self.head)
        self.assertEqual(newhead[1].split()[5], self.new[8])
        self.assertEqual(newhead[6].split()[3], self.new[0])
        self.assertEqual(newhead[6].split()[0], 'SIMU_CF3')

class SfacTableTest(unittest.TestCase):
    def setUp(self):
        self.res_file = './collect_resi.res'
        self.res_list = ResList(self.res_file)
        self.reslist =  self.res_list.get_res_list()
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

class NumberSchemeTest(unittest.TestCase):
    def setUp(self):
        self.numbers = ['O1A', 'C1A', 'C2A', 'F1A', 'F2A', 'F3A', 'C3A',
                        'F4A', 'F5A', 'F6A', 'C4A', 'F7A', 'F8A', 'F9A']
        res_file = './p21c.res'
        invert = True
        resi = False
        rl = ResList(res_file)
        reslist = rl.get_res_list()
        gdb = global_DB(invert)
        fragment = 'OC(cf3)3'
        dbatoms = gdb.get_atoms_from_fragment(fragment)
        self.num = NumberScheme(reslist, dbatoms, resi)

    def testrun_get_numberscheme(self):
        numberscheme = self.num.get_fragment_number_scheme()
        self.assertListEqual(numberscheme, self.numbers)

class insertAfixTest(unittest.TestCase):
    def setUp(self):
        import db
        self.maxDiff = None
        self.res_file = 'p21c.res'
        testresfile = './p21c.res'
        invert = True
        self.rl = ResList(testresfile)
        self.reslist = self.rl.get_res_list()
        self.dsrp = DSR_Parser(self.reslist, self.rl)
        self.dsr_dict = self.dsrp.parse_dsr_line()
        self.find_atoms = FindAtoms(self.reslist)
        self.gdb = global_DB(invert)
        fragment = 'OC(CF3)3'
        self.dbatoms = self.gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
        self.dbhead = self.gdb.get_head_from_fragment(fragment)        # this is only executed once
        self.resi = True #gdb.get_resi_from_fragment(fragment)
        self.dbtypes = get_atomtypes(self.dbatoms)
        with open('./intern.TXT') as txt:
            self.intern = txt.read()
        with open('./extern.TXT') as txt2:
            self.extern = txt2.read()
        misc.remove_file('dsr_TEST_p21c.res')
        self.sf = SfacTable(self.reslist, self.dbtypes)
        self.sfac_table = self.sf.set_sfac_table()
        self.num = NumberScheme(self.reslist, self.dbatoms, self.resi)
        self.numberscheme = self.num.get_fragment_number_scheme()
        self.db_testhead = db.db_testhead


    def testrun_afix(self):
        self.maxDiff = None
        afix = InsertAfix(self.reslist, self.dbatoms, self.dbtypes, self.dbhead, \
                          self.dsr_dict, self.sfac_table, self.find_atoms, self.numberscheme)
        afix_extern_entry = afix.build_afix_entry(True, self.res_file, 'TEST')
        afix_intern_entry = afix.build_afix_entry(False, self.res_file, 'TEST')
        self.assertEqual(afix_intern_entry, self.intern)
        self.assertEqual(afix_extern_entry, self.extern)
        misc.remove_file('dsr_TEST_p21c.res')

class removeDublicatesAfixTest(unittest.TestCase):
    def setUp(self):
        #self.verbosity = 4
        self.res_file = './collect_resi.res'
        self.res_list = ResList(self.res_file)
        self.reslist =  self.res_list.get_res_list()
        self.find_atoms = FindAtoms(self.reslist)
        invert = False
        self.dsrp = DSR_Parser(self.reslist, self.res_list)
        self.dsr_dict = self.dsrp.parse_dsr_line()
        fragment = 'OC(cf3)3'
        self.gdb = global_DB(invert)
        self.dbatoms = self.gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
        self.dbhead = self.gdb.get_head_from_fragment(fragment)        # this is only executed once
        self.dbtypes = get_atomtypes(self.dbatoms)
        #self.sf = SfacTable(self.reslist, self.dbtypes)
        #self.sfac_table = self.sf.set_sfac_table()
        self.sfac_table = ['C', 'H', 'N', 'O', 'F']
        self.resi = 'CCF3' #gdb.get_resi_from_fragment(fragment)
        self.num = NumberScheme(self.reslist, self.dbatoms, self.resi)
        self.numberscheme = self.num.get_fragment_number_scheme()
        self.afix = InsertAfix(self.reslist, self.dbatoms, self.dbtypes, self.dbhead, \
                          self.dsr_dict, self.sfac_table, self.find_atoms, self.numberscheme)
        self.db_testhead = ['SADI_CCF3 C1 C2 C1 C3 C1 C4', 'SADI_CCF3 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4 ',
                            'REM test']
        self.db_testhead2 = ['SADI_CF3 C1 C2 C1 C3 C1 C4 ', 'SADI_CF3 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4 ']

    def testrun_remove_dublicate_restraints(self):
        newhead = self.afix.remove_duplicate_restraints(self.db_testhead)#, self.resi.get_resiclass)
        newhead2 = self.afix.remove_duplicate_restraints(self.db_testhead2)#, self.resi.get_resiclass)
        self.assertListEqual(['', '', 'REM test'], newhead)
        self.assertListEqual(self.db_testhead2, newhead2)

    def testrun_remove_all_restraints(self):
        devided = [['SADI_CCF3 C1 C2 C1 C3 C1 C4\n', 'SADI_CCF3 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4\n'],
                   ['REM test\n']]
        newhead = self.afix.remove_all_restraints(self.db_testhead)
        self.assertListEqual(devided, newhead)

class invert_fragmentTest(unittest.TestCase):
    def setUp(self):
        self.dbatoms = [['O1', 3, '-0.01453', '1.66590', '0.10966'], ['C1', 1, '-0.00146', '0.26814', '0.06351'],
                        ['C2', 1, '-1.13341', '-0.23247', '-0.90730'], ['F1', 4, '-2.34661', '-0.11273', '-0.34544']]
        self.dbatoms2 = [['O1', 3, '-0.01453', '1.66590', '0.10966'], ['C1', 1, '-0.00146', '0.26814', '0.06351'],
                        ['C2', 1, '-1.13341', '-0.23247', '-0.90730'], ['F1', 4, '-2.34661', '-0.11273', '-0.34544']]
        self.inv_dbatoms = [['O1', 3, '0.01453', '-1.6659', '-0.10966'], ['C1', 1, '0.00146', '-0.26814', '-0.06351'],
                            ['C2', 1, '1.13341', '0.23247', '0.9073'], ['F1', 4, '2.34661', '0.11273', '0.34544']]

    def testrun_invert_dbatoms_coordinates(self):
        inverted = invert_dbatoms_coordinates(self.dbatoms)
        self.assertListEqual(self.inv_dbatoms, inverted)
        self.assertNotEqual(self.dbatoms2[0][2], inverted[0][2])

class atomsTest(unittest.TestCase):
    def setUp(self):
        self.atoms = Element.atoms

    def testrun_num_of_atoms(self):
        num_of_atoms = len(self.atoms)
        self.assertEqual(num_of_atoms, 92)
        self.assertNotEqual(num_of_atoms, 42)


class dbfileTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self._db_file_names = ("db1.txt", "db2.txt")
        self.rdb = ReadDB(dbdir='c:/test/')
        self.testnames = ['c:/test/db1.txt', 'c:/test/db2.txt']
        self.klein = ['\n', '<DMX>\n', 'REM test\n', 'RESI 3 TST1\n', 'SIMU C1\n','FRAG 17 1 1 1 90 90 90\n',
                      'O1  1  -1.3542148   -0.4780990   -0.5279749\n', '</DMX>']

    def testrun_dbpath(self):
        names = []
        for name in self._db_file_names:
            names.append(self.rdb.getDBpath(name))
        self.assertListEqual(names, self.testnames)

    def testrun_db_files_dict(self):
        db_file_names = ["db1_klein.TXT", "db2_klein.TXT"]
        rdb = ReadDB(dbdir='./', dbnames = db_file_names)
        self.assertListEqual(rdb.getDB_files_dict()['db2_klein'], self.klein)

    def testrun_find_db_tags(self):
        result = [['DME', 1, 'db1_klein'], ['DMX', 2, 'db2_klein']]
        db_file_names = ("db1_klein.TXT", "db2_klein.TXT")
        rdb = ReadDB(dbdir='./', dbnames = db_file_names)
        tags = rdb.find_db_tags()
        self.assertListEqual(result, tags)

    def testrun_dublicate_tags(self):
        db_file_names = ["db1_dublicate.TXT"]
        rdb = ReadDB(dbdir='./', dbnames = db_file_names)
        with self.assertRaises(SystemExit):
            rdb.find_db_tags()

    def testrun_ReadDB(self):
        db_file_names = ("db1_klein.TXT", "db2_klein.TXT")
        result = [['DME', 1, 'db1_klein'], ['DMX', 2, 'db2_klein']]
        rdb = ReadDB()
        rdb2 = ReadDB(dbdir='./', dbnames = db_file_names)
        names = rdb.find_db_tags()
        names2 = rdb2.find_db_tags()
        self.assertEqual(names, [])
        self.assertEqual(result, names2)

class globalDB(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.klein = ['RESI ClBE ', 'SADI C1 Cl1 C2 Cl2 ',
                      'SADI 0.04 C6 Cl1 C2 Cl1 C1 Cl2 C3 Cl2 ',
                      'SAME C2 > C6 C1 ', 'FLAT C1 > Cl2 ',
                      'SIMU C1 > Cl2 ', 'RIGU C1 > Cl2']
        self.result = {'dmx': {'comment': ['test'],
                               'head': ['SIMU C1'],
                               'resi': 'TST1',
                               'name': 'DMX',
                               'line': 2,
                               'db': 'db2_klein',
                               'fragline': ['FRAG', '17', '1', '1', '1', '90', '90', '90'],
                               'atoms': [['O1', '1', '-1.3542148', '-0.4780990', '-0.5279749']]},
                       'dme': {'comment': ['test'],
                               'head': [''],  # why is this an empty string if resi is False?
                               'resi': False,
                               'name': 'DME',
                               'line': 1,
                               'db': 'db1_klein',
                               'fragline': ['FRAG', '17', '1', '1', '1', '90', '90', '90'],
                               'atoms': [['O1', '1', '-1.3542148', '-0.4780990', '-0.5279749']]}}

    def testrun_get_head_lines(self):
        db_file_names = ("db1.TXT", "db2_klein.TXT")
        getdb = ReadDB(dbdir='.', dbnames = db_file_names)
        db_tags = getdb.find_db_tags()
        db_plain_dict = getdb.getDB_files_dict()
        #print(db_plain_dict)
        self.assertLessEqual(db_tags[0], ['DME-FREE', 1, 'db1'])
        gdb = global_DB(dbdir='.', dbnames = db_file_names)
        headlines = gdb.get_head_lines('DMe-free', 'db1', '1')
        #print(headlines)
        head = (['RESI DME', 'DFIX 1.409 O1 C1 O2 C4', 'DFIX 1.412 O1 C2 O2 C3',
          'DFIX 1.510 C2 C3', 'DANG 2.354 C1 C2 C3 C4', 'DANG 2.390 C2 O2 O1 C3',
          'SIMU O1 > C4', 'RIGU O1 > C4'],
         'FRAG 17 1 1 1 90 90 90',
         [['Name:', '1,2-Dimethoxyethane,', 'not', 'coordinated,', 'C4H10O2,', 'DME'],
          ['Src:', 'Turbomole,', 'B3-LYP/def2-TZVPP'],
          ['This', 'DME', 'is', 'not', 'coordinated']])
        self.assertListEqual(head[0], ['RESI DME', 'DFIX 1.409 O1 C1 O2 C4', 'DFIX 1.412 O1 C2 O2 C3',
          'DFIX 1.510 C2 C3', 'DANG 2.354 C1 C2 C3 C4', 'DANG 2.390 C2 O2 O1 C3',
          'SIMU O1 > C4', 'RIGU O1 > C4'])
        self.assertEqual(head[1], 'FRAG 17 1 1 1 90 90 90')
        self.assertListEqual(head[2], [['Name:', '1,2-Dimethoxyethane,', 'not', 'coordinated,', 'C4H10O2,', 'DME'],
          ['Src:', 'Turbomole,', 'B3-LYP/def2-TZVPP'],
          ['This', 'DME', 'is', 'not', 'coordinated']])


    def testrun_build_db_dict(self):
        db_file_names = ("db1_klein.TXT", "db2_klein.TXT")
        gdb = global_DB(dbdir='.', dbnames = db_file_names)
        db = gdb.build_db_dict()
        self.assertEqual(db['dmx']['line'], 2)
        self.assertEqual(db['dmx']['name'], 'DMX')
        self.assertEqual(db['dmx']['db'], 'db2_klein')
        self.assertEqual(db['dmx']['head'], ['SIMU C1'])
        self.assertEqual(db['dmx']['fragline'], ['FRAG', '17', '1', '1', '1', '90', '90', '90'])
        self.assertEqual(db['dmx']['resi'], 'TST1')
        self.assertEqual(db, self.result)

    def testrun_get_residue_from_head(self):
        db_file_names = ("db1_klein.TXT", "db2_klein.TXT")
        gdb = global_DB(dbdir='./', dbnames = db_file_names)
        gdb.build_db_dict()
        self.assertEqual(gdb.get_residue_from_head(self.klein), 'CLBE' )

    def testrun_get_residue_from_head2(self):
        # raises System exit, because residue in db_resinum.TXT is badly defined.
        db_file_names = ("db1_klein.TXT", "db2_klein.TXT", 'db_resinum.TXT')
        with self.assertRaises(SystemExit):
            gdb = global_DB(dbdir='./', dbnames = db_file_names)
            gdb.build_db_dict()

    def testrun_get_fragment_atoms(self):
        x = '-1.154'
        z = '0.526'
        o1 = ['O1', '1', '-1.154', '-0.748', '0.526']
        db_file_names = ("db1.TXT", "db2.TXT")
        gdb = global_DB(dbdir='./', dbnames = db_file_names)
        gdb.build_db_dict()
        atom = gdb.get_fragment_atoms('dme', 'db2', 1)[0]
        self.assertListEqual(o1, atom)
        self.assertEqual(x, atom[2])
        self.assertEqual(z, atom[4])
        self.assertEqual('1', atom[1])
        self.assertEqual('O1', atom[0])

    def testrun_get_fragment_atoms_shortline(self):
        db_file_names = ["db1_shortline.TXT"]
        gdb = global_DB(dbdir='./', dbnames = db_file_names)
        #db = gdb.build_db_dict()
        atom = gdb.get_fragment_atoms('dme-free', 'db1_shortline', 1)
        self.assertEqual(len(atom), 5)
        self.assertEqual(atom[0][0], 'O2')


    def testrun_get_fragment_atoms_inv(self):
        x = '1.154'
        z = '-0.526'
        o1 = ['O1', '1', '1.154', '0.748', '-0.526']
        db_file_names = ("db1.TXT", "db2.TXT")
        gdb = global_DB(invert = True, dbdir='./', dbnames = db_file_names)
        gdb.build_db_dict()
        atom = gdb.get_fragment_atoms('dme', 'db2', 1)[0]
        self.assertListEqual(o1, atom)
        self.assertEqual(x, atom[2])
        self.assertEqual(z, atom[4])
        self.assertEqual('1', atom[1])
        self.assertEqual('O1', atom[0])


    def testrun_get_fragment_atoms_noatoms(self):
        db_file_names = ("db1_noatoms.TXT", "db2.TXT")
        with self.assertRaises(SystemExit):
            gdb = global_DB(dbdir='./', dbnames = db_file_names)
            gdb.build_db_dict()
            gdb.get_fragment_atoms('dme-free', 'db1_noatoms', 1)


    def testrun_get_fragment_atoms_noend(self):
        db_file_names = ("db1_noend.TXT", "db2.TXT")
        with self.assertRaises(SystemExit):
            gdb = global_DB(dbdir='./', dbnames = db_file_names)
            gdb.build_db_dict()
            gdb.get_fragment_atoms('dme-free', 'db1_noend', 1)


    def testrun_header_consistency(self):
        self.maxDiff = None
        db_file_names = ["db1_head_inconsistent.TXT"]
        with self.assertRaises(SystemExit):
            gdb = global_DB(invert = True, dbdir='./', dbnames = db_file_names)
            db = gdb.build_db_dict()
            fragment = 'dmel'
            head = db[fragment]['head']
            gdb.check_db_header_consistency(head, fragment)

    def testrun_header_consistency2(self):
        self.maxDiff = None
        db_file_names = ["db1_head_inconsistent2.TXT"]
        with self.assertRaises(SystemExit):
            gdb = global_DB(invert = True, dbdir='./', dbnames = db_file_names)
            db = gdb.build_db_dict()
            fragment = 'dmem'
            head = db[fragment]['head']
            gdb.check_db_header_consistency(head, fragment)


    def testrun_get_comment_from_fragment1(self):
        self.maxDiff = None
        db_file_names = ["comment.TXT"]
        gdb = global_DB(invert = True, dbdir='./', dbnames = db_file_names)
        gdb.build_db_dict()
        fragment = 'com'
        comment = gdb.get_comment_from_fragment('com4')
        self.assertEqual(comment, 'A really fancy name.')
        names = ['name!', 'Name 1,2-Dimethoxyethane, not coordinated, C4H10O2, '
                 'DME, Src: Turbomole, B3-LYP/def2-TZVPP, This DME is not coordinated',
                 'Src: Turbomole, B3-LYP/def2-TZVPP, blub, This DME is not coordinated',
                 'A really fancy name.']
        for i in range(1,5):
            com = gdb.get_comment_from_fragment(fragment+str(i))
            self.assertEqual(com, names[i-1])



    def testrun_get_resi_from_fragment(self):
        self.maxDiff = None
        db_file_names = ["comment.TXT"]
        gdb = global_DB(invert = True, dbdir='./', dbnames = db_file_names)
        gdb.build_db_dict()
        fragment = 'com1'
        resi = gdb.get_resi_from_fragment(fragment)
        line = gdb.get_line_number_from_fragment(fragment)
        self.assertEqual(resi, 'DME')
        self.assertEqual(line, 2)


class ImportGRADE_Test(unittest.TestCase):
    def setUp(self):
        self.ig = ImportGRADE('./test-data/PFA.gradeserver_all.tgz')
        self.igi = ImportGRADE('./test-data/PFA.gradeserver_all.tgz', invert=True)

    # test for PFA1 is already in db and we want to import again

    def testrun_get_gradefiles(self):
        '''
        files[0] = pdb
        files[1] = dfix
        files[2] = obprop
        '''
        self.maxDiff = None
        files = self.ig.get_gradefiles()
        filenames = ['./grade-PFA.pdb', './grade-PFA.dfix', './obprop.txt']
        endings = []
        for num, i in enumerate(filenames):
            with open(i) as test_file:
                endings.append(test_file.readlines())
        self.assertListEqual(endings[num], files[num])

    def testrun_get_name_from_obprop(self):
        self.maxDiff = None
        filename = './obprop.txt'
        with open(filename) as filen:
            ob = filen.readlines()
        name = self.ig.get_name_from_obprop(ob)
        self.assertEqual(name, 'PFA')
        filename2 = './obprop_2.txt'
        with open(filename2) as filen:
            ob = filen.readlines()
        name = self.ig.get_name_from_obprop(ob)
        self.assertEqual(name, 'NONE')
        filename3 = './obprop_3.txt'
        with open(filename3) as filen:
            ob = filen.readlines()
        name = self.ig.get_name_from_obprop(ob)
        self.assertEqual(name, 'NONE')

    def testrun_get_comments(self):
        self.maxDiff = None
        filename = './grade-comments.dfix'
        ob = []
        with open(filename) as filen:
            for line in filen:
                #line = str(line, encoding='ascii')
                ob.append(line.split())
        comments = self.ig.get_comments()
        name = u'REM Name: ' + u'PFA'
        ob.insert(0, name.split()) # Name is always at first place
        self.assertEqual(comments, ob)

    def testrun_get_firstlast(self):
        files = self.ig.get_gradefiles()
        atoms = self.ig.get_pdbatoms(files[0])
        fl = self.ig.get_first_last_atom(atoms)
        self.assertTupleEqual(fl, ('AL1', 'F36'))

    def testrun_deleted_pdb_file(self):
        with self.assertRaises(SystemExit):
            ImportGRADE('./PFA.gradeserver_all_2.tgz')

    def testrun_get_restaraints(self):
        self.maxDiff = None
        restr = self.ig.get_restraints()
        # to generate the test file:
        #with open('test.txt', 'wb+') as file:
        #    for line in restr:
        #        file.write(' '.join(line)+'\n')
        filename = './grade_restraints.txt'
        tst = []
        with open(filename) as test_file:
            for line in test_file:
                tst.append(line.split())
        self.assertListEqual(restr, tst)

    def testrun_get_pdbatoms(self):
        #pdblines = []
        with open('./grade-PFA.pdb') as pdb_file:
            pdblines = pdb_file.readlines()
            pdbatoms = self.ig.get_pdbatoms(pdblines)
            self.assertListEqual(['AL1', '9.463', '-3.351', '3.397'], pdbatoms[0])

    def testrun_bild_grade_db_entry(self):
        import db
        dbentry = self.ig.bild_grade_db_entry()
        self.maxDiff = None
        dbtest = db.dbtest
        self.assertDictEqual(dbentry, dbtest)

class DSRParseTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.res_file = 'dsrparse.res'
        testresfile = './dsrparse.res'
        self.rl = ResList(testresfile)
        self.reslist = self.rl.get_res_list()
        self.dsrp = DSR_Parser(self.reslist, self.rl)
        self.dsr_dict = self.dsrp.parse_dsr_line()

    def testrun_find_dsr_command(self):
        self.maxDiff = None
        num = self.dsrp.find_dsr_command(line=False)
        self.assertEqual(num, 264)

class DSRParse2Test(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.res_file = 'dsrparse.res'
        testresfile = './dsrparse.res'
        self.rl = ResList(testresfile)
        self.reslist = self.rl.get_res_list()
        self.dsrp = DSR_Parser(self.reslist, self.rl)
        self.dsr_dict = self.dsrp.parse_dsr_line()

    def testrun_find_dsr_command(self):
        self.maxDiff = None
        line = self.dsrp.find_dsr_command(line=True)
        #self.assertEqual(num, 264)
        string = 'rem dsr put oc(cf3)3 with o1 c1 c2 c3 c4 on O1_3 c1_3 q6 Q4 q7 resi cf3  PART 2 occ -31 dfix\n'
        self.assertEqual(string, line)

class ExportTest(unittest.TestCase):
    def setUp(self):
        self.invert = False
        self.gdb = global_DB(self.invert)
        self.export_clip = 'benzene'

    def testrun_format_calced_coords(self):
        export = Export(self.export_clip, self.gdb)
        bigcell = export.format_calced_coords(['1', '1', '1', '90', '90', '90'])
        smallcell = export.format_calced_coords(['2', '1', '1', '90', '90', '90'])
        cell1 = ['50', '50', '50', '90', '90', '90']
        cell2 = ['2', '1', '1', '90', '90', '90']
        self.assertListEqual(bigcell, cell1)
        self.assertListEqual(smallcell, cell2)


    def testrun_export_to_clip(self):
        '''
        Exports the current fragment to the clipboard.
        '''
        from export import Export
        gdb = global_DB(self.invert)
        export = Export(self.export_clip, gdb)
#        with self.assertRaises(SystemExit):
        self.assertTrue(export.export_to_clip())

    def testrun_do_export_fragment(self):
        self.maxDiff = None
        fragment = 'toluene'
        from export import Export
        export = Export(fragment, self.gdb, self.invert)
        resfile = export.export_resfile()
        resgood = ['TITL toluene\n', 'REM This file was exported by DSR version {}\n'.format(VERSION),
                   'REM Name: Toluene, C7H8\nREM Source: GRADE import\nREM Gradeserver from http://grade.globalphasing.org\n',
                   'CELL 0.71073    50.000   50.000   50.000   90.000   90.000   90.000\n',
                   'ZERR    1.00   0.000    0.000    0.000    0.000    0.000    0.000\n', 'LATT  -1\n', 'SFAC C\n',
                   'UNIT 1 \n', 'WGHT  0.1\n', 'FVAR  1\n', 'rem Restraints from DSR database:\n',
                   'DFIX 1.509 0.011 C1 C2\nDFIX 1.385 0.012 C2 C3\nDFIX 1.385 0.012 C2 C7\nDFIX 1.384 0.010 C3 C4\nDFIX 1.379 0.015 C4 C5\nDFIX 1.379 0.015 C5 C6\nDFIX 1.384 0.010 C6 C7\nDANG 2.520 0.017 C1 C3\nDANG 2.520 0.017 C1 C7\nDANG 2.374 0.018 C3 C7\nDANG 2.410 0.019 C2 C4\nDANG 2.396 0.019 C3 C5\nDANG 2.386 0.022 C4 C6\nDANG 2.396 0.019 C5 C7\nDANG 2.410 0.019 C2 C6\nFLAT C2 C1 C3 C7\nFLAT C2 C3 C4 C5\nFLAT C3 C4 C5 C6\nFLAT C4 C5 C6 C7\nFLAT C5 C6 C7 C2\nFLAT C6 C7 C2 C3\nFLAT C7 C2 C3 C4\n',
                   '\n\n', 'C1    1     0.024000    -0.000460     0.072300   11.00   0.04\nC2    1     0.024060    -0.000240     0.042120   11.00   0.04\nC3    1     0.000300    -0.000220     0.027800   11.00   0.04\nC4    1     0.000300    -0.000020     0.000100   11.00   0.04\nC5    1     0.024160     0.000160    -0.013760   11.00   0.04\nC6    1     0.047960     0.000120     0.000180   11.00   0.04\nC7    1     0.047880    -0.000080     0.027880   11.00   0.04',
                   '\nHKLF 4\nEND\n']
        self.assertListEqual(resfile, resgood)

    def testrun_do_export_fragment_all(self):
        self.maxDiff = None
        fragment = 'toluene'
        from export import Export
        export = Export(fragment, self.gdb, self.invert, export_all=True)
        resfile = export.export_resfile()
        resgood = ['TITL toluene\n', 'REM This file was exported by DSR version {}\n'.format(VERSION),
                   'REM Name: Toluene, C7H8\nREM Source: GRADE import\nREM Gradeserver from http://grade.globalphasing.org\n',
                   'CELL 0.71073    50.000   50.000   50.000   90.000   90.000   90.000\n',
                   'ZERR    1.00   0.000    0.000    0.000    0.000    0.000    0.000\n', 'LATT  -1\n', 'SFAC C\n',
                   'UNIT 1 \n', 'WGHT  0.1\n', 'FVAR  1\n',
                   '\n\n', 'C1    1     0.024000    -0.000460     0.072300   11.00   0.04\nC2    1     0.024060    -0.000240     0.042120   11.00   0.04\nC3    1     0.000300    -0.000220     0.027800   11.00   0.04\nC4    1     0.000300    -0.000020     0.000100   11.00   0.04\nC5    1     0.024160     0.000160    -0.013760   11.00   0.04\nC6    1     0.047960     0.000120     0.000180   11.00   0.04\nC7    1     0.047880    -0.000080     0.027880   11.00   0.04',
                   '\nHKLF 4\nEND\n']
        self.assertListEqual(resfile, resgood)


class ResidueTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.res_file = './p21c.res'
        self.rl = ResList(self.res_file)
        self.res_list = self.rl.get_res_list()
        self.find_atoms = FindAtoms(self.res_list)
        self.rle = ResListEdit(self.res_list, self.find_atoms)
        self.dsrp = DSR_Parser(self.res_list, self.rle)
        self.dsr_dict = self.dsrp.parse_dsr_line()
        #fragment = self.dsr_dict['fragment']
        fragment = 'ch2cl2'
        invert = False
        self.head2 = ['SIMU C1 > CL2',
                 'RIGU C1 > CL2',
                 'DFIX 1.7602 0.002 CL1 C1 CL2 C1',
                 'DFIX 2.9338 0.004 CL1 CL2']
        self.gdb = global_DB(invert)
        self.residue_class = self.gdb.get_resi_from_fragment(fragment)
        self.db = self.gdb.build_db_dict()
        self.fragline = self.gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
        self.dbatoms = self.gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
        self.dbhead = self.gdb.get_head_from_fragment(fragment)
        self.resi = Resi(self.res_list, self.dsr_dict, self.dbhead, self.residue_class, self.find_atoms)


    def testrun_get_resinumber(self):
        self.assertEqual(self.resi.get_resinumber, '4')
        self.assertNotEqual(self.resi.get_resinumber, 'False')
        self.assertEqual(self.resi.get_residue_class, 'CF3')
        self.assertNotEqual(self.resi.get_residue_class, 'CCF3')


    def testrun_remove_resi(self):
        head = self.resi.remove_resi(self.dbhead)
        self.assertEqual(head, self.head2)


    def testrun_format_restraints(self):
        resihead1 = self.resi.format_restraints(self.head2)
        resihead2 = ['SIMU_CF3 C1 > CL2', 'RIGU_CF3 C1 > CL2',
                    'DFIX_CF3 1.7602 0.002 CL1 C1 CL2 C1',
                    'DFIX_CF3 2.9338 0.004 CL1 CL2']
        self.assertEqual(resihead1, resihead2)

    def testrun_get_unique_resinumber(self):
        num = self.resi.get_unique_resinumber('2')
        self.assertEqual(num, '4')
        self.assertNotEqual(num, '2')


    def testrun_get_resi_syntax(self):
        empty_dict = {'alias': None, 'class': None, 'number': None}
        only_number = {'alias': None, 'class': None, 'number': '2'}
        class_number = {'alias': None, 'class': 'CF3', 'number': '2'}
        only_class = {'alias': None, 'class': 'CF3', 'number': None}
        class_number_alias = {'alias': '3', 'class': 'CF3', 'number': '2'}
        self.assertDictEqual(empty_dict, self.resi.get_resi_syntax('RESI'))
        with self.assertRaises(SystemExit):
            self.assertDictEqual(empty_dict, self.resi.get_resi_syntax([]))
        with self.assertRaises(SystemExit):
            self.resi.get_resi_syntax(['dsfg', 'rgasr'])
        with self.assertRaises(SystemExit):
            self.resi.get_resi_syntax(['dgh3', '3', '23435'])
        with self.assertRaises(SystemExit):
            self.resi.get_resi_syntax(['dgh3', '23435'])
        self.assertDictEqual(self.resi.get_resi_syntax('2'.split()), only_number)
        self.assertEqual(self.resi.get_resi_syntax(['2', 'CF3']), class_number)
        self.assertEqual(self.resi.get_resi_syntax(['CF3','2', '3']), class_number_alias)
        self.assertEqual(self.resi.get_resi_syntax(['CF3']), only_class)


    def testrun_build_up_residue(self):
        dsr_dict1 = {'target': ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7'], 'fragment': 'OC(CF3)3',
                     'occupancy': '-31', 'source': ['O1', 'C1', 'C2', 'C3', 'C4'],
                     'resi': ['CF13', '2'], 'command': 'PUT', 'dfix': False, 'part': '2'}
        dsr_dict2 = {'target': ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7'], 'fragment': 'OC(CF3)3',
                     'occupancy': '-31', 'source': ['O1', 'C1', 'C2', 'C3', 'C4'],
                     'resi': ['CF23', '5'], 'command': 'PUT', 'dfix': False, 'part': '2'}
        resi1 = Resi(self.res_list, dsr_dict1, self.dbhead, self.residue_class, self.find_atoms)
        resi2 = Resi(self.res_list, dsr_dict2, self.dbhead, self.residue_class, self.find_atoms)
        residue1 = {'alias': None, 'class': 'CF13', 'number': '4'}
        residue2 = {'alias': None, 'class': 'CF23', 'number': '5'}
        self.assertDictEqual(resi1.build_up_residue(), residue1)
        self.assertDictEqual(resi2.build_up_residue(), residue2)


    def testrun_make_resihead(self):
        testhead = ['RESI TST', 'SADI 0.02 C1 C2 C1 C3 C1 C4', 'SADI 0.02 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4',
                   'SADI 0.04 C2 C3 C3 C4 C2 C4', 'SIMU O1 > F9', 'RIGU O1 > F9']
        dsr_dict1 = {'target': ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7'],
                     'fragment': 'OC(CF3)3',
                     'occupancy': '-31',
                     'source': ['O1', 'C1', 'C2', 'C3', 'C4'],
                     'resi': ['CF13'], # only class given
                     'command': 'PUT',
                     'dfix': False,
                     'part': '2'}
        dsr_dict2 = {'target': ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7'],
                     'fragment': 'OC(CF3)3',
                     'occupancy': '-31',
                     'source': ['O1', 'C1', 'C2', 'C3', 'C4'],
                     'resi': False, # no resi given
                     'command': 'PUT',
                     'dfix': False,
                     'part': '2'}
        dsr_dict3 = {'target': ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7'],
                     'fragment': 'OC(CF3)3',
                     'occupancy': '-31',
                     'source': ['O1', 'C1', 'C2', 'C3', 'C4'],
                     'resi': '', # only resi active, class from db
                     'command': 'PUT',
                     'dfix': False,
                     'part': '2'}
        dsr_dict4 = {'target': ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7'],
                     'fragment': 'OC(CF3)3',
                     'occupancy': '-31',
                     'source': ['O1', 'C1', 'C2', 'C3', 'C4'],
                     'resi': ['tes1', '8'], # class and number given
                     'command': 'PUT',
                     'dfix': False,
                     'part': '2'}
        dsr_dict5 = {'target': ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7'],
                     'fragment': 'OC(CF3)3',
                     'occupancy': '-31',
                     'source': ['O1', 'C1', 'C2', 'C3', 'C4'],
                     'resi': ['89'], # class and number given
                     'command': 'PUT',
                     'dfix': False,
                     'part': '2'}
        dsr_dict6 = {'target': ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7'],
                     'fragment': 'OC(CF3)3',
                     'occupancy': '-31',
                     'source': ['O1', 'C1', 'C2', 'C3', 'C4'],
                     'resi': ['tes1', '8', '10'], # class, number and alias given
                     'command': 'PUT',
                     'dfix': False,
                     'part': '2'}
        resi1 = Resi(self.res_list, dsr_dict1, testhead, 'CF13', self.find_atoms)
        resi2 = Resi(self.res_list, dsr_dict2, testhead, 'TEST', self.find_atoms)
        resi3 = Resi(self.res_list, dsr_dict3, testhead, 'TEST', self.find_atoms)
        resi4 = Resi(self.res_list, dsr_dict4, testhead, 'TES1', self.find_atoms)
        resi5 = Resi(self.res_list, dsr_dict5, testhead, 'TES1', self.find_atoms)
        resi6 = Resi(self.res_list, dsr_dict6, testhead, 'TES1', self.find_atoms)
        head1 = ['SADI_CF13 0.02 C1 C2 C1 C3 C1 C4', 'SADI_CF13 0.02 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4', 'SADI_CF13 0.04 C2 C3 C3 C4 C2 C4',
                'SIMU_CF13 O1 > F9', 'RIGU_CF13 O1 > F9', 'RESI 4 CF13']
        head2 = ['SADI 0.02 C1 C2 C1 C3 C1 C4', 'SADI 0.02 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4', 'SADI 0.04 C2 C3 C3 C4 C2 C4', 'SIMU O1 > F9', 'RIGU O1 > F9']
        head3 = ['SADI_TEST 0.02 C1 C2 C1 C3 C1 C4', 'SADI_TEST 0.02 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4', 'SADI_TEST 0.04 C2 C3 C3 C4 C2 C4',
                'SIMU_TEST O1 > F9', 'RIGU_TEST O1 > F9', 'RESI 4 TEST']
        head4 = ['SADI_TES1 0.02 C1 C2 C1 C3 C1 C4', 'SADI_TES1 0.02 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4', 'SADI_TES1 0.04 C2 C3 C3 C4 C2 C4',
                'SIMU_TES1 O1 > F9', 'RIGU_TES1 O1 > F9', 'RESI 8 TES1']
        head5 = ['SADI_TES1 0.02 C1 C2 C1 C3 C1 C4', 'SADI_TES1 0.02 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4', 'SADI_TES1 0.04 C2 C3 C3 C4 C2 C4',
                'SIMU_TES1 O1 > F9', 'RIGU_TES1 O1 > F9', 'RESI 89 TES1']
        head6 = ['SADI_TES1 0.02 C1 C2 C1 C3 C1 C4', 'SADI_TES1 0.02 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4', 'SADI_TES1 0.04 C2 C3 C3 C4 C2 C4',
                'SIMU_TES1 O1 > F9', 'RIGU_TES1 O1 > F9', 'RESI 10 TES1 8']
        # dsr command with residue
        self.assertListEqual(resi1.make_resihead(), head1)
        # dsr command without residue
        self.assertListEqual(resi2.make_resihead(), head2)
        # residue from database
        self.assertListEqual(resi3.make_resihead(), head3)
        # class and number given
        self.assertListEqual(resi4.make_resihead(), head4)
        # only number given
        self.assertListEqual(resi5.make_resihead(), head5)
        # class, number alias given
        self.assertListEqual(resi6.make_resihead(), head6)


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
        noatoms = misc.get_atoms([])
        self.assertEqual(noatoms, [])
        atoms = misc.get_atoms(self.dbatoms)
        self.assertListEqual(atoms, self.dbatoms)
        stratoms = misc.get_atoms(self.strdbatoms)
        self.assertListEqual(stratoms, self.dbatoms)
        tst = misc.get_atoms(['O1 3 -0.01453 1.66590 0.10966'])
        tst2 = misc.get_atoms(['O1 3 -0.01453 1.66590'])
        tst3 = misc.get_atoms('O1 3 -0.01453 1.66590  0.10966')
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
                    'SADI 0.1 F1 C1 F2 C1 F3 C1 F4 C1 F5 C1 F6 C1 F7 C1 F8 C1 F9 C1',
                    'SIMU O1 > F9', 'RIGU O1 > F9']
        head = misc.wrap_headlines(db_testhead)
        self.assertListEqual(head, wraphead)
        unwrap = misc.unwrap_head_lines(wraphead)
        self.assertListEqual(unwraped, unwrap)


    def testrun_makelist(self):
        lst = misc.makelist('Hallo Daniel!')
        self.assertListEqual(lst, ['HALLO', 'DANIEL!'])
        self.assertNotEqual(lst, ['HALLO', 'DANIEL'])


    def testrun_which(self):
        which = misc.which('twunk_32')[0]
        self.assertEqual(which, 'C:\\Windows\\twunk_32.exe')


    def testrun_zero(self):
        zero = misc.zero(3, 3)
        zer0 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        self.assertListEqual(zero, zer0)


    def testrun_matrixmult(self):
        m1 = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        m2 = [[2.01, 0, 0], [0, 2, 0], [0, 0, 2]]
        m3 = [[1.000001, 1.0, 1.0], [1.0, 2, 1.0], [1.01, 1.00, 2]]
        m4 = [[4.02, 0, 0], [0, 4, 0], [0.0, 0, 4]]
        m5 = [[2.000002, 2.0, 2.0], [2.0, 4.0, 2.0], [2.02, 2.0, 4.0]]
        mult1 = misc.matrix_mult(m1, m2)
        self.assertListEqual(mult1, m4)
        mult2 = misc.matrix_mult(m1, m3)
        self.assertListEqual(mult2, m5)


    def testrun_format_atom_names(self):
        names = misc.format_atom_names(atomnames, part=2, resinum=4)
        names2 = misc.format_atom_names(atomnames, part='2', resinum='4')
        names3 = misc.format_atom_names(atomnames, part='2')
        names4 = misc.format_atom_names(atomnames, resinum=91)
        names5 = misc.format_atom_names(atomnames, part='2a')
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
        cellf = [ float(i) for i in cell ]
        dst = misc.atomic_distance(coord1, coord2, cellf)
        self.assertAlmostEqual(6.05249787959, dst, 5)


    def testrun_frac_to_cart(self):
        coord1 = (-0.186843,   0.282708,   0.526803)
        cellf = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
        cart = (-2.741505423999065, 5.909586678000002, 10.775200700893734)
        N1 = misc.frac_to_cart(coord1, cellf)
        self.assertAlmostEqual(cart, N1, 8)

    def testrun_determinante(self):
        m1 = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        det = misc.determinante(m1)
        self.assertEqual(det, 8)

    def testrun_subtract(self):
        m1 = [2, 0, 0]
        m2 = [1, 1.5, 0.104]
        sub = misc.subtract_vect(m1, m2)
        self.assertEqual(sub, (1, -1.5, -0.104))


    def testrun_vol_tetrahedron(self):
        #cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
        a = (0.838817,   0.484526,   0.190081) # a ist um 0.01 ausgelenkt
        b = (0.875251,   0.478410,   0.256955)
        c = (0.789290,   0.456520,   0.301616)
        d = (0.674054,   0.430194,   0.280727)
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
        self.assertEqual(dice1, 0.837838)
        self.assertEqual(dice2, 0.789474)
        self.assertEqual(dice3, 1.0)
        self.assertEqual(dice4, 1.0)




if __name__ == "__main__":
    unittest.main()
