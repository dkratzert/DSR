#/usr/bin/env python
#-*- encoding: utf-8 -*-


import unittest
from atomhandling import get_atomtypes, FindAtoms, check_source_target,\
    rename_dbhead_atoms, SfacTable, Elem_2_Sfac, NumberScheme
from resfile import ResList
from dbfile import global_DB, invert_dbatoms_coordinates, ReadDB
from afix import InsertAfix
from dsrparse import DSR_Parser
import misc
from resi import Resi
from atoms import Element
import os



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
        self.res_file = 'unit-tests/collect_resi.res'
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
        self.assertEqual(self.fa.get_atomcoordinates(atoms_fe).keys()[0], fe_upper)
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
                    print(i.strip('\n'))
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
        self.assertEquals(newhead[1].split()[5], self.new[8])
        self.assertEquals(newhead[6].split()[3], self.new[0])
        self.assertEquals(newhead[6].split()[0], 'SIMU_CF3')

class SfacTableTest(unittest.TestCase):
    def setUp(self):
        self.res_file = 'unit-tests/collect_resi.res'
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
        res_file = 'unit-tests/p21c.res'
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
        self.res_file = 'p21c.res'
        testresfile = './unit-tests/p21c.res'
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
        with open('unit-tests/intern.TXT') as txt:
            self.intern = txt.read()
        with open('unit-tests/extern.TXT') as txt2:
            self.extern = txt2.read()
        misc.remove_file('TEST_p21c.res')
        self.sf = SfacTable(self.reslist, self.dbtypes)
        self.sfac_table = self.sf.set_sfac_table()
        self.num = NumberScheme(self.reslist, self.dbatoms, self.resi)
        self.numberscheme = self.num.get_fragment_number_scheme()
        self.db_testhead = ['SADI C1 C2 C1 C3 C1 C4 ', 'SADI F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4 ',
                       'SADI 0.04 C2 C3 C3 C4 C2 C4', 'SADI 0.04 O1 C2 O1 C3 O1 C4 ',
                       'SADI 0.04 F1 F2 F2 F3 F3 F1 F4 F5 F5 F6 F6 F4 F7 F8 F8 F9 F9 F7 ',
                       'SADI 0.1 F1 C1 F2 C1 F3 C1 F4 C1 F5 C1 F6 C1 F7 C1 F8 C1 F9 C1 ',
                       'SIMU O1 > F9 ', 'RIGU O1 > F9']


    def testrun_afix(self):
        afix = InsertAfix(self.reslist, self.dbatoms, self.dbtypes, self.dbhead, \
                          self.dsr_dict, self.sfac_table, self.find_atoms, self.numberscheme)
        afix_extern_entry = afix.build_afix_entry(True, self.res_file, 'TEST')
        afix_intern_entry = afix.build_afix_entry(False, self.res_file, 'TEST')
        self.assertEqual(afix_intern_entry, self.intern)
        self.assertEqual(afix_extern_entry, self.extern)
        misc.remove_file('TEST_p21c.res')

class removeDublicatesAfixTest(unittest.TestCase):
    def setUp(self):
        #self.verbosity = 4
        self.res_file = 'unit-tests/collect_resi.res'
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
        self._db_file_names = ("db1.txt", "db2.txt")
        self.rdb = ReadDB(dbdir='c:/test/')
        self.testnames = ['c:/test/db1.txt', 'c:/test/db2.txt']
        self.klein = ('<DMX>\n', 'REM test\n', 'FRAG 17 1 1 1 90 90 90\n',
                      'O1  1  -1.3542148   -0.4780990   -0.5279749\n', '</DMX>')

    def testrun_dbpath(self):
        names = []
        for name in self._db_file_names:
            names.append(self.rdb.getDBpath(name))
        self.assertListEqual(names, self.testnames)

    def testrun_db_files_dict(self):
        db_file_names = ("db1_klein.TXT", "db2_klein.TXT")
        rdb = ReadDB(dbdir='./unit-tests', dbnames = db_file_names)
        self.assertTupleEqual(rdb.getDB_files_dict()['db2_klein'], self.klein)

    def testrun_find_db_tags(self):
        result = [['DME', '1', 'db1_klein'], ['DMX', '1', 'db2_klein']]
        db_file_names = ("db1_klein.TXT", "db2_klein.TXT")
        rdb = ReadDB(dbdir='./unit-tests', dbnames = db_file_names)
        tags = rdb.find_db_tags()
        self.assertListEqual(result, tags)

    def testrun_dublicate_tags(self):
        db_file_names = ["db1_dublicate.TXT"]
        rdb = ReadDB(dbdir='./unit-tests', dbnames = db_file_names)
        with self.assertRaises(SystemExit):
            rdb.find_db_tags()


class globalDB(unittest.TestCase):
    def setUp(self):
        self.result = {'dmx': {'comment': ['test'], 'head': [''], 'resi': False, 'name': 'DMX', 'line': '1', 'db': 'db2_klein', 'fragline': ['FRAG', '17', '1', '1', '1', '90', '90', '90'], 'atoms': [['O1', '1', '-1.3542148', '-0.4780990', '-0.5279749']]}, 'dme': {'comment': ['test'], 'head': [''], 'resi': False, 'name': 'DME', 'line': '1', 'db': 'db1_klein', 'fragline': ['FRAG', '17', '1', '1', '1', '90', '90', '90'], 'atoms': [['O1', '1', '-1.3542148', '-0.4780990', '-0.5279749']]}}

    def testrun_build_db_dict(self):
        db_file_names = ("db1_klein.TXT", "db2_klein.TXT")
        #self.rdb = ReadDB(dbdir='./unit-tests', dbnames = db_file_names)
        #self.dbnames = self.rdb.find_db_tags()
#        invert = True
        gdb = global_DB(dbdir='./unit-tests', dbnames = db_file_names)
        db = gdb.build_db_dict()
        self.assertEqual(db, self.result)








if __name__ == "__main__":
    unittest.main()
