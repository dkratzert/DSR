#/usr/bin/env python
#-*- encoding: utf-8 -*-


import unittest

from afix import InsertAfix
from atomhandling import get_atomtypes, FindAtoms, check_source_target, \
    rename_dbhead_atoms, SfacTable, Elem_2_Sfac, NumberScheme
from atoms import Element
from dbfile import global_DB, invert_dbatoms_coordinates, ReadDB, ImportGRADE
from dsrparse import DSR_Parser
import misc
from resfile import ResList


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
        self.maxDiff = None
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
        self.maxDiff = None
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
        rdb = ReadDB(dbdir='./unit-tests', dbnames = db_file_names)
        self.assertListEqual(rdb.getDB_files_dict()['db2_klein'], self.klein)

    def testrun_find_db_tags(self):
        result = [['DME', 1, 'db1_klein'], ['DMX', 2, 'db2_klein']]
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

    def testrun_build_db_dict(self):
        db_file_names = ("db1_klein.TXT", "db2_klein.TXT")
        #self.rdb = ReadDB(dbdir='./unit-tests', dbnames = db_file_names)
        #self.dbnames = self.rdb.find_db_tags()
#        invert = True
        gdb = global_DB(dbdir='./unit-tests', dbnames = db_file_names)
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
        gdb = global_DB(dbdir='./unit-tests', dbnames = db_file_names)
        gdb.build_db_dict()
        self.assertEqual(gdb.get_residue_from_head(self.klein), 'CLBE' )

    def testrun_get_residue_from_head2(self):
        db_file_names = ("db1_klein.TXT", "db2_klein.TXT", 'db_resinum.TXT')
        with self.assertRaises(SystemExit):
            gdb = global_DB(dbdir='./unit-tests', dbnames = db_file_names)
            gdb.build_db_dict()

    def testrun_get_fragment_atoms(self):
        x = '-1.154'
        z = '0.526'
        o1 = ['O1', '1', '-1.154', '-0.748', '0.526']
        db_file_names = ("db1.TXT", "db2.TXT")
        gdb = global_DB(dbdir='./unit-tests', dbnames = db_file_names)
        gdb.build_db_dict()
        atom = gdb.get_fragment_atoms('dme', 'db2', 1)[0]
        self.assertListEqual(o1, atom)
        self.assertEqual(x, atom[2])
        self.assertEqual(z, atom[4])
        self.assertEqual('1', atom[1])
        self.assertEqual('O1', atom[0])

    def testrun_get_fragment_atoms_shortline(self):
        db_file_names = ["db1_shortline.TXT"]
        gdb = global_DB(dbdir='./unit-tests', dbnames = db_file_names)
        #db = gdb.build_db_dict()
        atom = gdb.get_fragment_atoms('dme-free', 'db1_shortline', 1)
        self.assertEqual(len(atom), 5)
        self.assertEqual(atom[0][0], 'O2')


    def testrun_get_fragment_atoms_inv(self):
        x = '1.154'
        z = '-0.526'
        o1 = ['O1', '1', '1.154', '0.748', '-0.526']
        db_file_names = ("db1.TXT", "db2.TXT")
        gdb = global_DB(invert = True, dbdir='./unit-tests', dbnames = db_file_names)
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
            gdb = global_DB(dbdir='./unit-tests', dbnames = db_file_names)
            gdb.build_db_dict()
            gdb.get_fragment_atoms('dme-free', 'db1_noatoms', 1)


    def testrun_get_fragment_atoms_noend(self):
        db_file_names = ("db1_noend.TXT", "db2.TXT")
        with self.assertRaises(SystemExit):
            gdb = global_DB(dbdir='./unit-tests', dbnames = db_file_names)
            gdb.build_db_dict()
            gdb.get_fragment_atoms('dme-free', 'db1_noend', 1)


    def testrun_header_consistency(self):
        self.maxDiff = None
        db_file_names = ["db1_head_inconsistent.TXT"]
        with self.assertRaises(SystemExit):
            gdb = global_DB(invert = True, dbdir='./unit-tests', dbnames = db_file_names)
            db = gdb.build_db_dict()
            fragment = 'dmel'
            head = db[fragment]['head']
            gdb.check_db_header_consistency(head, fragment)

    def testrun_header_consistency2(self):
        self.maxDiff = None
        db_file_names = ["db1_head_inconsistent2.TXT"]
        with self.assertRaises(SystemExit):
            gdb = global_DB(invert = True, dbdir='./unit-tests', dbnames = db_file_names)
            db = gdb.build_db_dict()
            fragment = 'dmem'
            head = db[fragment]['head']
            gdb.check_db_header_consistency(head, fragment)


    def testrun_get_comment_from_fragment1(self):
        self.maxDiff = None
        db_file_names = ["comment.TXT"]
        gdb = global_DB(invert = True, dbdir='./unit-tests', dbnames = db_file_names)
        gdb.build_db_dict()
        fragment = 'com'
        names = ['name!', 'Name 1,2-Dimethoxyethane, not coordinated, C4H10O2, '
                 'DME, Src: Turbomole, B3-LYP/def2-TZVPP, This DME is not coordinated',
                 'Src: Turbomole, B3-LYP/def2-TZVPP, blub, This DME is not coordinated',
                 'A really fancy name.']
        for i in range(1,5):
            com = gdb.get_comment_from_fragment(fragment+str(i))
            #print(com)
            self.assertEqual(com, names[i-1])

    def testrun_get_resi_from_fragment(self):
        self.maxDiff = None
        db_file_names = ["comment.TXT"]
        gdb = global_DB(invert = True, dbdir='./unit-tests', dbnames = db_file_names)
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


    def testrun_get_gradefiles(self):
        '''
        files[0] = pdb
        files[1] = dfix
        files[2] = obprop
        '''
        self.maxDiff = None
        files = self.ig.get_gradefiles()
        filenames = ['unit-tests/grade-PFA.pdb', 'unit-tests/grade-PFA.dfix', 'unit-tests/obprop.txt']
        endings = []
        for num, i in enumerate(filenames):
            with open(i) as test_file:
                #file = file.decode('ascii')
                endings.append(test_file.readlines())
        self.assertListEqual(endings[num], files[num])

    def testrun_get_name_from_obprop(self):
        self.maxDiff = None
        filename = 'unit-tests/obprop.txt'
        with open(filename) as filen:
            ob = filen.readlines()
        name = self.ig.get_name_from_obprop(ob)
        self.assertEqual(name, 'PFA')
        filename2 = 'unit-tests/obprop_2.txt'
        with open(filename2) as filen:
            ob = filen.readlines()
        name = self.ig.get_name_from_obprop(ob)
        self.assertEqual(name, 'NONE')
        filename3 = 'unit-tests/obprop_3.txt'
        with open(filename3) as filen:
            ob = filen.readlines()
        name = self.ig.get_name_from_obprop(ob)
        self.assertEqual(name, 'NONE')

    def testrun_get_comments(self):
        self.maxDiff = None
        filename = 'unit-tests/grade-comments.dfix'
        ob = []
        with open(filename) as filen:
            for line in filen:
                ob.append(line.split())
        comments = self.ig.get_comments()
        name = 'REM Name: ' + 'PFA'
        ob.append(name.split())
        self.assertEqual(comments, ob)

    def testrun_get_firstlast(self):
        files = self.ig.get_gradefiles()
        atoms = self.ig.get_pdbatoms(files[0])
        fl = self.ig.get_first_last_atom(atoms)
        self.assertTupleEqual(fl, (u'AL1', u'F36'))

    def testrun_deleted_pdb_file(self):
        with self.assertRaises(SystemExit):
            ImportGRADE('./unit-tests/PFA.gradeserver_all_2.tgz')

    def testrun_get_restaraints(self):
        self.maxDiff = None
        restr = self.ig.get_restraints()
        # to generate the test file:
        #with open('test.txt', 'wb+') as file:
        #    for line in restr:
        #        file.write(' '.join(line)+'\n')
        filename = './unit-tests/grade_restraints.txt'
        tst = []
        with open(filename) as test_file:
            for line in test_file:
                tst.append(line.split())
        self.assertListEqual(restr, tst)

    def testrun_get_pdbatoms(self):
        #pdblines = []
        with open('unit-tests/grade-PFA.pdb') as pdb_file:
            pdblines = pdb_file.readlines()
            pdbatoms = self.ig.get_pdbatoms(pdblines)
            self.assertListEqual([u'AL1', u'9.463', u'-3.351', u'3.397'], pdbatoms[0])

    def testrun_bild_grade_db_entry(self):
        dbentry = self.ig.bild_grade_db_entry()
        dbtest = {'PFA1': {'comment': [[u'REM', u'Produced', u'by', u'Grade', u'Web', u'Server', u'http://grade.globalphasing.org'],
                                       [u'REM', u'GEN:', u'Generated', u'by', u'GRADE', u'1.2.7', u'(February', u'19', u'2014)'],
                                       [u'REM', u'GEN:', u'from', u'mol2', u'file'],
                                       [u'REM', u'GEN:', u'using', u'quantum', u'mechanics', u'PM3'],
                                       [u'REM', u'grade-cif2shelx', u'output'],
                                       [u'REM', u'grade-cif2shelx', u'extracts', u'restraints', u'from', u'a', u'grade', u'CIF', u'file'],
                                       [u'REM', u'Version:', u'0.0.5', u'<Dec', u'20', u'2013>'], [u'REM', u'Name:', u'PFA']],
                           'head': [['DFIX', '1.785', '0.030', 'AL1', 'O1'], ['DFIX', '1.784', '0.030', 'AL1', 'O2'],
                                    ['DFIX', '1.788', '0.030', 'AL1', 'O3'], ['DFIX', '1.785', '0.030', 'AL1', 'O4'],
                                    ['DFIX', '1.325', '0.030', 'C1', 'O1'], ['DFIX', '1.620', '0.030', 'C1', 'C2'],
                                    ['DFIX', '1.628', '0.030', 'C1', 'C3'], ['DFIX', '1.622', '0.030', 'C1', 'C4'],
                                    ['DFIX', '1.344', '0.030', 'C2', 'F1'], ['DFIX', '1.352', '0.030', 'C2', 'F2'],
                                    ['DFIX', '1.347', '0.030', 'C2', 'F3'], ['DFIX', '1.347', '0.030', 'C3', 'F4'],
                                    ['DFIX', '1.351', '0.030', 'C3', 'F5'], ['DFIX', '1.345', '0.030', 'C3', 'F6'],
                                    ['DFIX', '1.344', '0.030', 'C4', 'F7'], ['DFIX', '1.348', '0.030', 'C4', 'F8'],
                                    ['DFIX', '1.352', '0.030', 'C4', 'F9'], ['DFIX', '1.329', '0.030', 'C5', 'O2'],
                                    ['DFIX', '1.630', '0.030', 'C5', 'C6'], ['DFIX', '1.624', '0.030', 'C5', 'C7'],
                                    ['DFIX', '1.622', '0.030', 'C5', 'C8'], ['DFIX', '1.351', '0.030', 'C6', 'F10'],
                                    ['DFIX', '1.347', '0.030', 'C6', 'F11'], ['DFIX', '1.345', '0.030', 'C6', 'F12'],
                                    ['DFIX', '1.348', '0.030', 'C7', 'F13'], ['DFIX', '1.352', '0.030', 'C7', 'F14'],
                                    ['DFIX', '1.344', '0.030', 'C7', 'F15'], ['DFIX', '1.352', '0.030', 'C8', 'F16'],
                                    ['DFIX', '1.347', '0.030', 'C8', 'F17'], ['DFIX', '1.344', '0.030', 'C8', 'F18'],
                                    ['DFIX', '1.324', '0.030', 'C9', 'O3'], ['DFIX', '1.626', '0.030', 'C10', 'C9'],
                                    ['DFIX', '1.622', '0.030', 'C11', 'C9'], ['DFIX', '1.620', '0.030', 'C12', 'C9'],
                                    ['DFIX', '1.347', '0.030', 'C10', 'F19'], ['DFIX', '1.345', '0.030', 'C10', 'F20'],
                                    ['DFIX', '1.351', '0.030', 'C10', 'F21'], ['DFIX', '1.345', '0.030', 'C11', 'F22'],
                                    ['DFIX', '1.348', '0.030', 'C11', 'F23'], ['DFIX', '1.352', '0.030', 'C11', 'F24'],
                                    ['DFIX', '1.352', '0.030', 'C12', 'F25'], ['DFIX', '1.347', '0.030', 'C12', 'F26'],
                                    ['DFIX', '1.345', '0.030', 'C12', 'F27'], ['DFIX', '1.327', '0.030', 'C13', 'O4'],
                                    ['DFIX', '1.622', '0.030', 'C13', 'C14'], ['DFIX', '1.344', '0.030', 'C14', 'F28'],
                                    ['DFIX', '1.348', '0.030', 'C14', 'F29'], ['DFIX', '1.352', '0.030', 'C14', 'F30'],
                                    ['DFIX', '1.622', '0.030', 'C13', 'C15'], ['DFIX', '1.629', '0.030', 'C13', 'C16'],
                                    ['DFIX', '1.345', '0.030', 'C15', 'F31'], ['DFIX', '1.348', '0.030', 'C15', 'F32'],
                                    ['DFIX', '1.352', '0.030', 'C15', 'F33'], ['DFIX', '1.351', '0.030', 'C16', 'F34'],
                                    ['DFIX', '1.347', '0.030', 'C16', 'F35'], ['DFIX', '1.345', '0.030', 'C16', 'F36'],
                                    ['DANG', '2.981', '0.062', 'O1', 'O2'], ['DANG', '2.918', '0.064', 'O1', 'O3'],
                                    ['DANG', '2.856', '0.066', 'O2', 'O3'], ['DANG', '2.866', '0.065', 'O1', 'O4'],
                                    ['DANG', '2.948', '0.063', 'O2', 'O4'], ['DANG', '2.920', '0.064', 'O3', 'O4'],
                                    ['DANG', '3.069', '0.044', 'AL1', 'C1'], ['DANG', '2.502', '0.054', 'C2', 'O1'],
                                    ['DANG', '2.408', '0.056', 'C3', 'O1'], ['DANG', '2.577', '0.062', 'C2', 'C3'],
                                    ['DANG', '2.501', '0.054', 'C4', 'O1'], ['DANG', '2.577', '0.062', 'C2', 'C4'],
                                    ['DANG', '2.578', '0.062', 'C3', 'C4'], ['DANG', '2.506', '0.055', 'C1', 'F1'],
                                    ['DANG', '2.478', '0.055', 'C1', 'F2'], ['DANG', '2.127', '0.055', 'F1', 'F2'],
                                    ['DANG', '2.472', '0.055', 'C1', 'F3'], ['DANG', '2.159', '0.054', 'F1', 'F3'],
                                    ['DANG', '2.138', '0.055', 'F2', 'F3'], ['DANG', '2.474', '0.056', 'C1', 'F4'],
                                    ['DANG', '2.492', '0.055', 'C1', 'F5'], ['DANG', '2.145', '0.055', 'F4', 'F5'],
                                    ['DANG', '2.505', '0.055', 'C1', 'F6'], ['DANG', '2.158', '0.054', 'F4', 'F6'],
                                    ['DANG', '2.129', '0.055', 'F5', 'F6'], ['DANG', '2.503', '0.055', 'C1', 'F7'],
                                    ['DANG', '2.480', '0.055', 'C1', 'F8'], ['DANG', '2.160', '0.054', 'F7', 'F8'],
                                    ['DANG', '2.482', '0.055', 'C1', 'F9'], ['DANG', '2.126', '0.055', 'F7', 'F9'],
                                    ['DANG', '2.138', '0.055', 'F8', 'F9'], ['DANG', '3.019', '0.046', 'AL1', 'C5'],
                                    ['DANG', '2.391', '0.057', 'C6', 'O2'], ['DANG', '2.519', '0.054', 'C7', 'O2'],
                                    ['DANG', '2.582', '0.062', 'C6', 'C7'], ['DANG', '2.519', '0.054', 'C8', 'O2'],
                                    ['DANG', '2.575', '0.062', 'C6', 'C8'], ['DANG', '2.573', '0.062', 'C7', 'C8'],
                                    ['DANG', '2.496', '0.055', 'C5', 'F10'], ['DANG', '2.477', '0.056', 'C5', 'F11'],
                                    ['DANG', '2.140', '0.055', 'F10', 'F11'], ['DANG', '2.505', '0.055', 'C5', 'F12'],
                                    ['DANG', '2.127', '0.055', 'F10', 'F12'], ['DANG', '2.158', '0.054', 'F11', 'F12'],
                                    ['DANG', '2.482', '0.055', 'C5', 'F13'], ['DANG', '2.479', '0.055', 'C5', 'F14'],
                                    ['DANG', '2.138', '0.055', 'F13', 'F14'], ['DANG', '2.505', '0.055', 'C5', 'F15'],
                                    ['DANG', '2.164', '0.054', 'F13', 'F15'], ['DANG', '2.126', '0.055', 'F14', 'F15'],
                                    ['DANG', '2.479', '0.055', 'C5', 'F16'], ['DANG', '2.471', '0.055', 'C5', 'F17'],
                                    ['DANG', '2.138', '0.055', 'F16', 'F17'], ['DANG', '2.511', '0.054', 'C5', 'F18'],
                                    ['DANG', '2.123', '0.055', 'F16', 'F18'], ['DANG', '2.162', '0.054', 'F17', 'F18'],
                                    ['DANG', '3.071', '0.044', 'AL1', 'C9'], ['DANG', '2.411', '0.056', 'C10', 'O3'],
                                    ['DANG', '2.521', '0.054', 'C11', 'O3'], ['DANG', '2.579', '0.062', 'C10', 'C11'],
                                    ['DANG', '2.472', '0.055', 'C12', 'O3'], ['DANG', '2.582', '0.062', 'C10', 'C12'],
                                    ['DANG', '2.574', '0.062', 'C11', 'C12'], ['DANG', '2.471', '0.056', 'C9', 'F19'],
                                    ['DANG', '2.503', '0.055', 'C9', 'F20'], ['DANG', '2.156', '0.054', 'F19', 'F20'],
                                    ['DANG', '2.490', '0.055', 'C9', 'F21'], ['DANG', '2.148', '0.054', 'F19', 'F21'],
                                    ['DANG', '2.129', '0.055', 'F20', 'F21'], ['DANG', '2.507', '0.055', 'C9', 'F22'],
                                    ['DANG', '2.480', '0.055', 'C9', 'F23'], ['DANG', '2.159', '0.054', 'F22', 'F23'],
                                    ['DANG', '2.479', '0.055', 'C9', 'F24'], ['DANG', '2.125', '0.055', 'F22', 'F24'],
                                    ['DANG', '2.139', '0.055', 'F23', 'F24'], ['DANG', '2.483', '0.055', 'C9', 'F25'],
                                    ['DANG', '2.470', '0.055', 'C9', 'F26'], ['DANG', '2.138', '0.055', 'F25', 'F26'],
                                    ['DANG', '2.502', '0.055', 'C9', 'F27'], ['DANG', '2.130', '0.055', 'F25', 'F27'],
                                    ['DANG', '2.161', '0.054', 'F26', 'F27'], ['DANG', '3.034', '0.045', 'AL1', 'C13'],
                                    ['DANG', '2.506', '0.055', 'C13', 'F28'], ['DANG', '2.471', '0.055', 'C13', 'F29'],
                                    ['DANG', '2.157', '0.054', 'F28', 'F29'], ['DANG', '2.484', '0.055', 'C13', 'F30'],
                                    ['DANG', '2.129', '0.055', 'F28', 'F30'], ['DANG', '2.141', '0.055', 'F29', 'F30'],
                                    ['DANG', '2.469', '0.055', 'C14', 'O4'], ['DANG', '2.541', '0.053', 'C15', 'O4'],
                                    ['DANG', '2.579', '0.062', 'C14', 'C15'], ['DANG', '2.410', '0.056', 'C16', 'O4'],
                                    ['DANG', '2.576', '0.062', 'C14', 'C16'], ['DANG', '2.577', '0.062', 'C15', 'C16'],
                                    ['DANG', '2.507', '0.055', 'C13', 'F31'], ['DANG', '2.480', '0.055', 'C13', 'F32'],
                                    ['DANG', '2.163', '0.054', 'F31', 'F32'], ['DANG', '2.477', '0.055', 'C13', 'F33'],
                                    ['DANG', '2.125', '0.055', 'F31', 'F33'], ['DANG', '2.138', '0.055', 'F32', 'F33'],
                                    ['DANG', '2.493', '0.055', 'C13', 'F34'], ['DANG', '2.477', '0.056', 'C13', 'F35'],
                                    ['DANG', '2.143', '0.055', 'F34', 'F35'], ['DANG', '2.506', '0.055', 'C13', 'F36'],
                                    ['DANG', '2.127', '0.055', 'F34', 'F36'], ['DANG', '2.158', '0.054', 'F35', 'F36'],
                                    [u'SIMU', u'AL1', u'>', u'F36'], [u'RIGU', u'AL1', u'>', u'F36']],
                           'resi': 'PFA1',
                           'name': 'PFA1',
                           'line': None,
                           'db': 'dsr-user-db',
                           'fragline': ['FRAG', '17', '1', '1', '1', '90', '90', '90'],
                           'atoms': [[u'AL1', u'9.463', u'-3.351', u'3.397'], [u'O1', u'8.422', u'-2.079', u'4.093'],
                                     [u'C1', u'7.600', u'-1.044', u'4.188'], [u'C2', u'6.402', u'-1.006', u'3.088'],
                                     [u'F1', u'6.808', u'-0.868', u'1.814'], [u'F2', u'5.537', u'0.018', u'3.271'],
                                     [u'F3', u'5.634', u'-2.113', u'3.126'], [u'C3', u'6.864', u'-1.113', u'5.645'],
                                     [u'F4', u'7.745', u'-1.050', u'6.663'], [u'F5', u'5.994', u'-0.102', u'5.866'],
                                     [u'F6', u'6.133', u'-2.224', u'5.844'], [u'C4', u'8.321', u'0.415', u'4.119'],
                                     [u'F7', u'9.299', u'0.599', u'5.022'], [u'F8', u'8.862', u'0.674', u'2.912'],
                                     [u'F9', u'7.473', u'1.448', u'4.340'], [u'O2', u'10.476', u'-2.774', u'2.046'],
                                     [u'C5', u'11.284', u'-2.995', u'1.015'], [u'C6', u'12.221', u'-1.663', u'0.860'],
                                     [u'F10', u'13.116', u'-1.734', u'-0.152'], [u'F11', u'11.498', u'-0.555', u'0.606'],
                                     [u'F12', u'12.972', u'-1.388', u'1.941'], [u'C7', u'10.554', u'-3.196', u'-0.430'],
                                     [u'F13', u'9.842', u'-4.339', u'-0.502'], [u'F14', u'11.424', u'-3.271', u'-1.464'],
                                     [u'F15', u'9.723', u'-2.199', u'-0.779'], [u'C8', u'12.317', u'-4.246', u'1.170'],
                                     [u'F16', u'13.197', u'-4.348', u'0.147'], [u'F17', u'13.085', u'-4.139', u'2.272'],
                                     [u'F18', u'11.742', u'-5.460', u'1.217'], [u'O3', u'8.442', u'-4.660', u'2.731'],
                                     [u'C9', u'7.678', u'-5.739', u'2.656'], [u'C10', u'6.896', u'-5.722', u'1.224'],
                                     [u'F19', u'7.746', u'-5.796', u'0.181'], [u'F20', u'6.140', u'-4.628', u'1.018'],
                                     [u'F21', u'6.038', u'-6.754', u'1.061'], [u'C11', u'8.439', u'-7.177', u'2.732'],
                                     [u'F22', u'9.376', u'-7.369', u'1.786'], [u'F23', u'9.047', u'-7.386', u'3.917'],
                                     [u'F24', u'7.606', u'-8.232', u'2.581'], [u'C12', u'6.524', u'-5.782', u'3.801'],
                                     [u'F25', u'5.670', u'-6.824', u'3.672'], [u'F26', u'5.736', u'-4.689', u'3.777'],
                                     [u'F27', u'6.987', u'-5.893', u'5.059'], [u'O4', u'10.460', u'-3.990', u'4.732'],
                                     [u'C14', u'12.313', u'-2.660', u'5.662'], [u'C13', u'11.237', u'-3.873', u'5.802'],
                                     [u'F28', u'11.776', u'-1.428', u'5.645'], [u'F29', u'13.057', u'-2.764', u'4.542'],
                                     [u'F30', u'13.213', u'-2.613', u'6.672'], [u'C15', u'10.516', u'-3.668', u'7.248'],
                                     [u'F31', u'9.650', u'-4.639', u'7.590'], [u'F32', u'9.843', u'-2.503', u'7.335'],
                                     [u'F33', u'11.391', u'-3.632', u'8.281'], [u'C16', u'12.124', u'-5.241', u'5.935'],
                                     [u'F34', u'13.010', u'-5.223', u'6.957'], [u'F35', u'11.358', u'-6.327', u'6.154'],
                                     [u'F36', u'12.878', u'-5.519', u'4.856']]}}
        self.assertDictEqual(dbentry, dbtest)

class DSRParseTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.res_file = 'dsrparse.res'
        testresfile = './unit-tests/dsrparse.res'
        invert = True
        self.rl = ResList(testresfile)
        self.reslist = self.rl.get_res_list()
        self.dsrp = DSR_Parser(self.reslist, self.rl)
        self.dsr_dict = self.dsrp.parse_dsr_line()

    #def testrun_


if __name__ == "__main__":
    unittest.main()
