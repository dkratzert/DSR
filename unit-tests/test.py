# /usr/bin/env python
# -*- encoding: utf-8 -*-


import doctest
import sys
import unittest
from os import system
import os

import dbfile

sys.path.insert(1, os.path.join(sys.path[0], '..'))
import afix
import atomhandling
import dsrparse
import elements
import export
import misc
import networkx
from afix import Afix
from atomhandling import get_atomtypes, FindAtoms, check_source_target, \
    rename_restraints_atoms, SfacTable, Elem_2_Sfac, NumberScheme
from atoms import Element, atoms
from dbfile import invert_atomic_coordinates, ImportGRADE
import dsr
from dsr import VERSION
from dsrparse import DSRParser
from export import Export
from resfile import ResList, ResListEdit
from resi import Resi
from restraints import format_atom_names

print(sys.version)


class doctestsTest(unittest.TestCase):
    def testrun_doctest(self):
        for name in [dsr, afix, dsrparse, export, misc, elements, networkx,
                     atomhandling, networkx.classes.graph, dsrparse, dbfile]:
            failed, attempted = doctest.testmod(name)  # , verbose=True)
            if failed == 0:
                print('passed all {} tests in {}!'.format(attempted, name.__name__))
            else:
                msg = '!!!!!!!!!!!!!!!! {} of {} tests failed in {}  !!!!!!!!!!!!!!!!!!!!!!!!!!!'.format(failed,
                                                                                                         attempted,
                                                                                                         name.__name__)
                self.assertFalse(failed, msg)


# @unittest.skip("skipping dsr_complete_runs_Test ")
class dsr_complete_runs_Test(unittest.TestCase):
    def setUp(self):
        self.maxDiff = 20
        # remove this to view the results:
        # self.dsr = '/Applications/DSR/dsr'
        # self.dsr = 'D:\Programme\DSR\dsr'
        self.dsr = r'"{}"'.format(misc.which('dsr')[-1])
        print(self.dsr)

        # 1 -r resi cf3 part 2 occ -31
        # 2 -r resi cf3 part 2 occ -31 dfix
        # 3 -r part 2 occ -31
        # 4 -re resi cf3 part 2 occ -31
        # 5 -re resi cf3 part 2 occ -31 dfix
        # 6 -re part 2 occ -31
        # 7 rigid
        # 8 -s
        # 9 -e toliene
        # 10 -r -t
        # 11 -r replace
        # 12 -r resi replace

    def dsr_runtest(self, nummer=99, parameter='-r', external_file='', hkl=1,
                    limit_start=6, limit_end=-1, ending='res', remlines=None):
        """
        runs a test where the whole dsr is started with different input files
        and compares the result with optimal output
        :type nummer: int
        :type external_file: sting
        :param nummer: res file number
        :param parameter: dsr command line parameter
        :return: None
        """
        if remlines is None:
            remlines = []
        d = []
        c = []
        print('{} '.format(nummer) * 10, 'start:')
        misc.copy_file('test-data/beispiel/{}.hkl'.format(hkl), 'test-data/beispiel/{}a.hkl'.format(nummer))
        misc.copy_file('test-data/beispiel/{}.res'.format(nummer), 'test-data/beispiel/{}a.res'.format(nummer))
        system('{0} {1} ./test-data/beispiel/{2}a.res'.format(self.dsr, parameter, nummer))
        with open('./test-data/beispiel/{}a.{}'.format(nummer, ending)) as txt:
            a = txt.readlines()[limit_start:limit_end]
        with open('./test-data/beispiel/{}-erg.{}'.format(nummer, ending)) as txt2:
            b = txt2.readlines()[limit_start:limit_end]
        if external_file:
            with open('./test-data/beispiel/{}.dfix'.format(external_file)) as ext:
                c = ext.readlines()
            with open('./test-data/beispiel/{}-erg.dfix'.format(external_file)) as ext2:
                d = ext2.readlines()
        print('{} test:'.format(nummer))
        print("parameter:", parameter)
        for line in remlines:
            del a[line]
            del b[line]
        self.assertEqual(b, a)
        if external_file:
            self.assertEqual(d, c)
        print('{} '.format(nummer) * 10, "ende")
        misc.remove_file('./test-data/beispiel/{}a.hkl'.format(nummer))
        misc.remove_file('./test-data/beispiel/{}a.ins'.format(nummer))
        misc.remove_file('./test-data/beispiel/{}a.res'.format(nummer))
        misc.remove_file('./test-data/beispiel/{}a.lst'.format(nummer))
        misc.remove_file('./test-data/beispiel/{}a.fcf'.format(nummer))
        misc.remove_file('./test-data/beispiel/{}.fcf'.format(nummer))
        misc.remove_file('./test-data/beispiel/{}.dfix'.format(external_file))

    # @unittest.skip(" skipping1 ")
    def testrun_run1(self):
        """
        regular dsr run with
        resi cf3 PART 2 occ -31
        """
        self.maxDiff = 10
        self.dsr_runtest(1, '-r')

    # @unittest.skip(" skipping2 ")
    def testrun_run2(self):
        """
        regular -r dsr run with
        resi cf3 dfix =
            PART 2 occ -31
        """
        self.dsr_runtest(2, '-r')

    # @unittest.skip(" skipping3 ")
    def testrun_run3(self):
        """
        regular run with:
         occ -31 PART 2
        """
        self.dsr_runtest(3, '-r')

    # @unittest.skip(" skipping4 ")
    def testrun_run4(self):
        """
        external -re restraints with:
        resi cf3 =
            PART 2 occ -31
        """
        self.dsr_runtest(4, '-re', external_file='dsr_CCF3_6_4a')

    # @unittest.skip(" skipping5 ")
    def testrun_run5(self):
        """
        dsr -re  with:
        REM dsr put oc(CF3)3 with o1 c1 c2 C3 on O1_4 C1_4 Q6 Q7 resi Ccf3 part 2 occ =
            -31 dfix
        """
        self.dsr_runtest(5, '-re', external_file='dsr_CCF3_6_5a_dfx')

    # @unittest.skip(" skipping6 ")
    def testrun_run6(self):
        """
        -re   PART 2 occ -31
        """
        self.dsr_runtest(6, '-re', external_file='dsr_1_6a')

    # @unittest.skip(" skipping 7")
    def testrun_run7(self):
        """
        rigid

        """
        self.dsr_runtest(7, '-g -r')

    # @unittest.skip(" skipping 8")
    def testrun_8(self):
        """
        dsr -s tol

        """
        print('8 ' * 10, 'start')
        system("{} -s tol > search.txt".format(self.dsr))
        with open('search.txt') as txt:
            se = txt.readlines()
        misc.remove_file('./search.txt')
        good = False
        for line in se:
            if line.startswith(" toluene           | Toluene, C7H8"):
                good = True
        self.assertTrue(good, "Search text differs")
        print('8 ' * 10, 'ende')

    # @unittest.skip(" skipping 9")
    def testrun_9(self):
        """
        dsr -e tol

        """
        print('9 ' * 10, 'start')
        print("{} -e tolUene".format(self.dsr))
        system("{} -e tolUene".format(self.dsr))
        with open('toluene.res') as txt:
            ex = txt.readlines()
        with open('toluene-erg.res') as txt2:
            ex_erg = txt2.readlines()
        del ex[1]  # line with the version number
        del ex_erg[1]
        self.assertEqual(ex, ex_erg)
        misc.remove_file('toluene.res')
        print('9 ' * 10, 'ende')

    # @unittest.skip(" skipping 10")
    def testrun_run10(self):
        """
        invert fragment
        -r -t
        """
        # Line 74 has one digit difference in Windows and Mac:
        self.dsr_runtest(10, '-t -r', limit_end=300, remlines=[74])

    # @unittest.skip(" skipping 11")
    def testrun_run11(self):
        """
        regular dsr run with
        replace resi PART 0
        """
        self.dsr_runtest(11, '-r', remlines=[108])

    # @unittest.skip(" skipping 12")
    def testrun_run12(self):
        """
        regular dsr run without fit with
        resi cf3 PART 2 occ -31
        """
        self.dsr_runtest(12, ending='ins', parameter='-n -r', remlines=[100, 115])

    # @unittest.skip(" skipping 13")
    def testrun_run13(self):
        """
        rem dsr put oc(cf3)3 with o1 c1 c2 c3 c4 on O1_3 c1_3 q6 Q4 q7 resi cf3 =
            PART 2 occ 501
        """
        self.dsr_runtest(13, '-r', remlines=[331])

    # @unittest.skip(" skipping 14")
    def testrun_run14(self):
        """
        dsr -r p21c.res
        REM  dsr put oc(CF3)3 with o1 c1 c2 C3 on O1_4 C1_4 Q6 Q7 resi cf3 =
            PART 2
        """
        self.dsr_runtest(14, '-r')

    # @unittest.skip(" skipping 15")
    def testrun_run15(self):
        """
        dsr -r foo
        REM  dsr put oc(CF3)3 with o1 c1 c2 C3 on O1_4 C1_4 Q6 Q7
        """
        self.dsr_runtest(15, '-r')

    # @unittest.skip(" skipping 16")
    def testrun_run16(self):
        """
        dsr -r
        REM  dsr put oc(CF3)3 with o1 c1 c2 C3 on O1_4 C1_4 Q6 Q7 RESI
        """
        self.dsr_runtest(16, '-r')

    # @unittest.skip(" skipping 17")
    def testrun_run17(self):
        """
        dsr -r with:
        REM  dsr put oc(CF3)3 with o1 c1 c2 C3 on O1_4 C1_4 Q6 Q7 resi occ -31
        """
        self.dsr_runtest(17, '-r')

    # @unittest.skip(" skipping 18")
    def testrun_run18(self):
        """
        dsr -r foo
        REM  dsr put oc(CF3)3 with o1 c1 c2 C3 on O1_4 C1_4 Q6 Q7 occ -31
        """
        self.dsr_runtest(18, '-r')

    # @unittest.skip(" skipping 19")
    def testrun_run19(self):
        """
        occ -31 without FVAR

        """
        self.dsr_runtest(19, '-re', external_file='dsr_1_19a')

    # @unittest.skip(" skipping 20")
    def testrun_run20(self):
        """
        rem dsr put CF6 on C22 split
        (AFIX 120)

        """
        self.dsr_runtest(20, '-r', hkl=20)

    # @unittest.skip(" skipping 21")
    def testrun_run21(self):
        """
        rem dsr put CF3 on C22
        (AFIX 130)
        """
        self.dsr_runtest(21, '-r', hkl=20, remlines=[27])

    # @unittest.skip(" skipping 22")
    def testrun_run22(self):
        """
        rem dsr put CF9 on C1

        """
        self.dsr_runtest(22, '-r', hkl=20)

    # @unittest.skip(" skipping 23")
    def testrun_run23(self):
        """
        REM DSR PUT TOLUENE WITH C2 C3 C5 ON Q1 Q2 Q1 PART -1 OCC 10.5 RESI TOL
        1.0005, 0.5447, 0.5342, 0.9314, 0.5395, 0.5126, 0.9995, 0.4553, 0.4658
        """
        self.dsr_runtest(23, '-target 1.0005 0.5447 0.5342 0.9314 0.5395 0.5126 0.9995 0.4553 0.4658 -r')


db_testhead = ['SADI C1 C2 C1 C3 C1 C4',
               'SADI F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4',
               'SADI 0.04 C2 C3 C3 C4 C2 C4', 'SADI 0.04 O1 C2 O1 C3 O1 C4',
               'SADI 0.04 F1 F2 F2 F3 F3 F1 F4 F5 F5 F6 F6 F4 F7 F8 F8 F9 F9 F7',
               'SADI 0.1 F1 C1 F2 C1 F3 C1 F4 C1 F5 C1 F6 C1 F7 C1 F8 C1 F9 C1',
               'SIMU O1 > F9', 'RIGU O1 > F9']

wraphead = ['SADI C1 C2 C1 C3 C1 C4\n',
            'SADI F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4 F4 C3 F5 C3 F6 =\n   C3 F7 C4 F8 C4 F9 C4\n',
            'SADI 0.04 C2 C3 C3 C4 C2 C4\n', 'SADI 0.04 O1 C2 O1 C3 O1 C4\n',
            'SADI 0.04 F1 F2 F2 F3 F3 F1 F4 F5 F5 F6 F6 F4 F7 F8 F8 F9 F9 F7\n',
            'SADI 0.1 F1 C1 F2 C1 F3 C1 F4 C1 F5 C1 F6 C1 F7 C1 F8 C1 F9 C1\n',
            'SIMU O1 > F9\n', 'RIGU O1 > F9\n']

atomnames = ['O1', 'C1', 'C2', 'F1', 'F2', 'F3', 'C3', 'F4', 'F5', 'F6', 'C4', 'F7', 'F8', 'F9']
formated_atoms = ['O1_4B', 'C1_4B', 'C2_4B', 'F1_4B', 'F2_4B', 'F3_4B', 'C3_4B', 'F4_4B', 'F5_4B', 'F6_4B', 'C4_4B',
                  'F7_4B', 'F8_4B', 'F9_4B']
part_atoms = ['O1_B', 'C1_B', 'C2_B', 'F1_B', 'F2_B', 'F3_B', 'C3_B', 'F4_B', 'F5_B', 'F6_B', 'C4_B', 'F7_B', 'F8_B',
              'F9_B']
resi_atoms = ['O1_91', 'C1_91', 'C2_91', 'F1_91', 'F2_91', 'F3_91', 'C3_91', 'F4_91', 'F5_91', 'F6_91', 'C4_91',
              'F7_91', 'F8_91', 'F9_91']
resi_atoms2 = ['O1_4', 'C1_4', 'C2_4', 'F1_4', 'F2_4', 'F3_4', 'C3_4', 'F4_4', 'F5_4', 'F6_4', 'C4_4', 'F7_4', 'F8_4',
               'F9_4']

coord1 = (0.100759, 0.449528, 0.438703)  # F4
coord1a = ('0.100759', '0.449528', '0.438703')  # F4
coord2 = (-0.362398, 0.278516, 0.447770)  # F10 6.052A
cell = [10.5086, 20.9035, 20.5072, 90, 94.13, 90]
cells = ['10.5086', '20.9035', '20.5072', '90', '94.13', '90']


class ElementsTest(unittest.TestCase):
    def setUp(self):
        self.el = Element()

    def testrun_get_atomic_number(self):
        num = self.el.get_atomic_number('Fe')
        ele = self.el.get_element(26)
        self.assertEqual(num, 26)
        self.assertEqual(ele, 'Fe')


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


class reslistTest(unittest.TestCase):
    def setUp(self):
        # print ("setUp reslistTest executed!")
        self.atom1 = 'F10   4   -0.362398   0.278516   0.447770  11.00000   0.02302   0.04023 =\n'
        self.atom2 = '0.02897   0.00131  -0.01216   0.00374\n'
        self.res_file = 'p21c.res'
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


class collect_residuesTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.res_file = './collect_resi.res'
        self.res_list = ResList(self.res_file)
        self.reslist = self.res_list.get_res_list()
        self.fa = FindAtoms(self.reslist)

    def testrun_collect_residues(self):
        # print(self.fa.collect_residues())
        collected_resi = {'0': [['FE1', [0.100001, 0.200002, 0.300003], 19, None, '0', 'FE', '5']],
                          '1': [['O1', [0.584527, 0.749093, 0.406892], 34, 'CCF3', 1, 'F', '3'],
                                ['H2', [0.2, 0.3, 0.4], 36, 'CCF3', 1, 'H', '6'],
                                ['C1', [0.462797, 0.766414, 0.415951], 37, 'CCF3', 1, 'C', '1'],
                                ['H1', [0.2, 0.3, 0.4], 40, 'CCF3', 1, 'H', '6']]}
        resi_null = [['FE1', [0.100001, 0.200002, 0.300003], 19, None, '0', 'FE', '5']]
        resi_fe = 'FE1'
        self.assertEqual(self.fa.collect_residues(), collected_resi)
        self.assertEqual(self.fa.collect_residues()['0'], resi_null)
        self.assertEqual(self.fa.collect_residues()['0'][0][0], resi_fe)

    def testrun_get_atomcoordinates(self):
        atoms_fe = ['Fe1']
        atoms_fe_0 = ['Fe1_0']
        atoms = ['Fe1', 'C1_1', 'O1_1']
        fe_upper = 'FE1'
        fe = {'FE1': [0.100001, 0.200002, 0.300003]}
        fe_str = {'FE1': ['0.100001', '0.200002', '0.300003']}
        coords = [0.100001, 0.200002, 0.300003]
        atoms_coords = {'C1_1': [0.462797, 0.766414, 0.415951],
                        'O1_1': [0.584527, 0.749093, 0.406892],
                        'FE1': [0.100001, 0.200002, 0.300003]}

        self.assertEqual(self.fa.get_atomcoordinates(atoms_fe), fe)
        self.assertNotEqual(self.fa.get_atomcoordinates(atoms_fe), fe_str)
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


class rename_DBHeadatomsTest(unittest.TestCase):
    def setUp(self):
        self.head = ['SADI_CF3 C1 C2 C1 C3 C1 C4', 'SADI_CF3 F1 C2 F2 C2 F3 C2 F4 C3'
                                                   ' F5 C3 F6 C3 F7 C4 F8 C4 F9 C4', 'SADI_CF3 0.04 C2 C3 C3 C4 C2 C4',
                     'SADI_CF3 0.04 O1 C2 O1 C3 O1 C4', 'SADI_CF3 0.04 F1 F2 F2 F3 F3'
                                                        ' F1 F4 F5 F5 F6 F6 F4 F7 F8 F8 F9 F9 F7',
                     'SADI_CF3 0.1 F1 C1 F2 '
                     'C1 F3 C1 F4 C1 F5 C1 F6 C1 F7 C1 F8 C1 F9 C1', 'SIMU_CF3 O1 > F9',
                     'RIGU_CF3 O1 > F9', 'RESI 4 CF3']
        self.old = ('O1', 'C1', 'C2', 'F1', 'F2', 'F3', 'C3', 'F4', 'F5', 'F6', 'C4', 'F7', 'F8', 'F9')
        self.new = ('F9A', 'F8A', 'F7A', 'C4A', 'F6A', 'F5A', 'F4A', 'C3A', 'F3A', 'F2A', 'F1A', 'C2A', 'C1A', 'O1A')
        self.new_rev = tuple(reversed(
            ['F9A', 'F8A', 'F7A', 'C4A', 'F6A', 'F5A', 'F4A', 'C3A', 'F3A', 'F2A', 'F1A', 'C2A', 'C1A', 'O1A']))

    def testrun_rename_dbheadatoms(self):
        newhead = rename_restraints_atoms(self.new, self.old, self.head)
        self.assertEqual(newhead[1].split()[5], self.new[8])
        self.assertEqual(newhead[6].split()[3], self.new[0])
        self.assertEqual(newhead[6].split()[0], 'SIMU_CF3')


class SfacTableTest(unittest.TestCase):
    def setUp(self):
        self.res_file = './collect_resi.res'
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


class NumberSchemeTest(unittest.TestCase):
    def setUp(self):
        self.numbers = ['O1A', 'C1A', 'C2A', 'F1A', 'F2A', 'F3A', 'C3A',
                        'F4A', 'F5A', 'F6A', 'C4A', 'F7A', 'F8A', 'F9A']
        res_file = './p21c.res'

        class Dsrp():
            resiflag = False

        dsrp = Dsrp()
        rl = ResList(res_file)
        reslist = rl.get_res_list()
        gdb = dbfile.ParseDB('../dsr_db.txt')
        fragment = 'OC(cf3)3'
        dbatoms = gdb.get_atoms(fragment)
        self.num = NumberScheme(reslist, dbatoms, dsrp)

    def testrun_get_numberscheme(self):
        numberscheme = self.num.get_fragment_number_scheme()
        self.assertListEqual(numberscheme, self.numbers)


@unittest.skip("skipping insertAfixTest")
class insertAfixTest(unittest.TestCase):
    def setUp(self):
        from options import OptionsParser
        import db
        self.maxDiff = None
        self.res_file = 'p21c.res'
        testresfile = './p21c.res'
        invert = True

        class OptionsParser():
            rigid_group = False
            target_coords = False

        self.options = OptionsParser()
        self.rl = ResList(testresfile)
        self.reslist = self.rl.get_res_list()
        self.dsrp = DSRParser(self.reslist)
        self.find_atoms = FindAtoms(self.reslist)
        self.gdb = dbfile.ParseDB()
        fragment = 'OC(CF3)3'
        self.dbatoms = self.gdb.get_atoms(fragment)  # only the atoms of the dbentry as list
        self.dbhead = self.gdb.get_restraints(fragment)  # this is only executed once
        self.resi = Resi(self.dsrp, 'CF3', self.find_atoms)
        self.dbtypes = get_atomtypes(self.dbatoms)
        with open('./intern.TXT') as txt:
            self.intern = txt.read()
        with open('./extern.TXT') as txt2:
            self.extern = txt2.read()
        misc.remove_file('./dsr_CF3_4_dsr_CF3_p21c.dfix')
        self.sf = SfacTable(self.reslist, self.dbtypes)
        self.sfac_table = self.sf.set_sfac_table()
        self.num = NumberScheme(self.reslist, self.dbatoms, self.dsrp)
        self.numberscheme = self.num.get_fragment_number_scheme()
        self.db_testhead = db.db_testhead


class removeDublicatesAfixTest(unittest.TestCase):
    def setUp(self):
        # self.verbosity = 4
        self.res_file = './collect_resi.res'
        self.res_list = ResList(self.res_file)
        self.reslist = self.res_list.get_res_list()
        self.find_atoms = FindAtoms(self.reslist)
        invert = False
        self.dsrp = DSRParser(self.reslist)
        fragment = 'OC(cf3)3'
        self.gdb = dbfile.ParseDB('../dsr_db.txt')
        self.dbatoms = self.gdb.get_atoms(fragment)  # only the atoms of the dbentry as list
        self.dbhead = self.gdb.get_restraints(fragment)  # this is only executed once
        self.dbtypes = get_atomtypes(self.dbatoms)
        # self.sf = SfacTable(self.reslist, self.dbtypes)
        # self.sfac_table = self.sf.set_sfac_table()
        self.sfac_table = ['C', 'H', 'N', 'O', 'F']
        self.resi = 'CCF3'  # gdb.get_resi_from_fragment(fragment)
        self.num = NumberScheme(self.reslist, self.dbatoms, self.dsrp)
        self.numberscheme = self.num.get_fragment_number_scheme()
        self.afix = Afix(self.reslist, self.dbatoms, self.dbtypes, self.dbhead,
                         self.dsrp, self.sfac_table, self.find_atoms, self.numberscheme, {'rigid_group': False})
        self.db_testhead = ['SADI_CCF3 C1 C2 C1 C3 C1 C4',
                            'SADI_CCF3 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4 ',
                            'REM test']
        self.db_testhead2 = ['SADI_CF3 C1 C2 C1 C3 C1 C4 ',
                             'SADI_CF3 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4 ']

    def testrun_remove_dublicate_restraints(self):
        newhead = self.afix.remove_duplicate_restraints(self.db_testhead,
                                                        self.afix.collect_all_restraints())  # , self.resi.get_resiclass)
        newhead2 = self.afix.remove_duplicate_restraints(self.db_testhead2,
                                                         self.afix.collect_all_restraints())  # , self.resi.get_resiclass)
        self.assertListEqual(['', '', 'REM test'], newhead)
        self.assertListEqual(self.db_testhead2, newhead2)


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
    def setUp(self):
        self.atoms = atoms

    def testrun_num_of_atoms(self):
        num_of_atoms = len(self.atoms)
        self.assertEqual(num_of_atoms, 92)
        self.assertNotEqual(num_of_atoms, 42)


class dbfileTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.rdb = dbfile.ParseDB('./db1.txt', "db2.txt")
        self.klein = ['\n', '<DMX>\n', 'REM test\n', 'RESI 3 TST1\n', 'SIMU C1\n', 'FRAG 17 1 1 1 90 90 90\n',
                      'O1  1  -1.3542148   -0.4780990   -0.5279749\n', '</DMX>']

    def testrun_dublicate_tags(self):
        with self.assertRaises(SystemExit):
            rdb = dbfile.ParseDB("db1_dublicate.TXT")


class globalDB(unittest.TestCase):
    def setUp(self):
        self.maxDiff = 9999
        self.klein = ['RESI ClBE ', 'SADI C1 Cl1 C2 Cl2 ',
                      'SADI 0.04 C6 Cl1 C2 Cl1 C1 Cl2 C3 Cl2 ',
                      'SAME C2 > C6 C1 ', 'FLAT C1 > Cl2 ',
                      'SIMU C1 > Cl2 ', 'RIGU C1 > Cl2']
        self.result = {'dmx': {'startline': 2,
                               'name': '',
                               'comments': ['REM test'],
                               'atoms': [['O1', 1, -1.3542148, -0.478099, -0.5279749]],
                               'cell': [1.0, 1.0, 1.0, 90.0, 90.0, 90.0],
                               'source': '',
                               'resi': 'TST1',
                               'restraints': ['SIMU C1'],
                               'endline': 8,
                               'dbname': 'dsr_user_db'},
                       'dme': {'startline': 1,
                               'name': '',
                               'comments': ['REM test'],
                               'atoms': [['O1', 1, -1.3542148, -0.478099, -0.5279749]],
                               'cell': [1.0, 1.0, 1.0, 90.0, 90.0, 90.0],
                               'source': '',
                               'resi': '',
                               'restraints': [],
                               'endline': 5,
                               'dbname': 'dsr_db'}
                       }
    
    def testrun_build_db_dict(self):
        gdb = dbfile.ParseDB("db1_klein.TXT", "db2_klein.TXT")
        db = gdb.databases
        self.assertEqual(db['dmx']['startline'], 2)
        self.assertEqual(db['dmx']['name'], '')
        self.assertEqual(db['dmx']['dbname'], 'dsr_user_db')
        self.assertEqual(db['dmx']['restraints'], ['SIMU C1'])
        self.assertEqual(db['dmx']['cell'], [1, 1, 1, 90, 90, 90])
        self.assertEqual(db['dmx']['resi'], 'TST1')
        self.assertEqual(db, self.result)

    def testrun_get_residue_from_head2(self):
        # raises System exit, because residue in db_resinum.TXT is badly defined.
        main_dbpath = "db_resinum.TXT"
        user_dbpath = "db1.TXT"
        with self.assertRaises(SystemExit):
            gdb = dbfile.ParseDB(main_dbpath, user_dbpath)

    def testrun_get_fragment_atoms(self):
        x = -1.154
        z = 0.526
        o1 = ['O1', 1, -1.154, -0.748, 0.526]
        gdb = dbfile.ParseDB("db1.TXT", "db2.TXT")
        atom = gdb.get_atoms('dme')[0]
        self.assertListEqual(o1, atom)
        self.assertEqual(x, atom[2])
        self.assertEqual(z, atom[4])
        self.assertEqual(1, atom[1])
        self.assertEqual('O1', atom[0])

    def testrun_get_fragment_atoms_shortline(self):
        gdb = dbfile.ParseDB("db1_shortline.TXT")
        # db = gdb.build_db_dict()
        atom = gdb.get_atoms('dme-free')
        self.assertEqual(len(atom), 5)
        self.assertEqual(atom[0][0], 'O2')

    def testrun_get_fragment_atoms_inv(self):
        x = -1.154
        z = 0.526
        o1 = ['O1', 1, -1.154, -0.748, 0.526]
        gdb = dbfile.ParseDB("db1.TXT", "db2.TXT")
        atom = gdb.get_atoms('dme')[0]
        self.assertListEqual(o1, atom)
        self.assertEqual(x, atom[2])
        self.assertEqual(z, atom[4])
        self.assertEqual(1, atom[1])
        self.assertEqual('O1', atom[0])

    def testrun_get_fragment_atoms_noatoms(self):
        with self.assertRaises(SystemExit):
            gdb = dbfile.ParseDB("db1_noatoms.TXT", "db2.TXT")
            gdb.get_atoms('dme-free')

    def testrun_get_fragment_atoms_noend(self):
        with self.assertRaises(SystemExit):
            gdb = dbfile.ParseDB("db1_noend.TXT", "db2.TXT")
            gdb.get_atoms('dme-free')

    def testrun_header_consistency(self):
        self.maxDiff = None
        main_dbpath = "./db1_head_inconsistent.TXT"
        user_dbpath = "./db2_klein.TXT"
        with self.assertRaises(SystemExit):
            gdb = dbfile.ParseDB(main_dbpath, user_dbpath)
            fragment = 'dmel'
            gdb.check_db_header_consistency(fragment)

    def testrun_header_consistency2(self):
        self.maxDiff = None
        maindb = "./db1_head_inconsistent2.TXT"
        userdb = "./db2_klein.TXT"
        with self.assertRaises(SystemExit):
            gdb = dbfile.ParseDB(maindb, userdb)
            fragment = 'dmem'
            head = gdb[fragment]['restraints']
            atoms = gdb.get_atoms(fragment, invert=True)
            gdb.check_db_header_consistency(fragment)

    def testrun_get_resi_from_fragment(self):
        self.maxDiff = None
        gdb = dbfile.ParseDB("comment.TXT", 'db1.txt')
        fragment = 'com1'
        resi = gdb.get_resi(fragment)
        line = gdb.get_startline(fragment)
        self.assertEqual(resi, 'DME')
        self.assertEqual(line, 2)


class ImportGRADE_Test(unittest.TestCase):
    def setUp(self):
        gdb = dbfile.ParseDB('../userdb.txt')
        self.ig = ImportGRADE('./test-data/PFA.gradeserver_all.tgz', gdb)
        self.igi = ImportGRADE('./test-data/PFA.gradeserver_all.tgz', gdb, invert=True)

    # test for PFA1 is already in db and we want to import again

    def testrun_get_gradefiles(self):
        '''
        files[0] = pdb
        files[1] = dfix
        files[2] = obprop
        '''
        self.maxDiff = None
        files = self.ig.get_gradefiles('./test-data/PFA.gradeserver_all.tgz')
        filenames = ['./grade-PFA.pdb', './grade-PFA.dfix']
        endings = []
        for num, i in enumerate(filenames):
            with open(i) as test_file:
                #   endings.append(test_file.readlines())
                tmp = []
                for line in test_file:
                    if not isinstance(line, str):
                        line = line.decode()
                    tmp.append(line)
                endings.append(tmp)
        self.assertListEqual(endings[num], files[num])

    def testrun_get_comments(self):
        self.maxDiff = None
        filename = './grade-comments.dfix'
        ob = []
        with open(filename) as filen:
            for line in filen:
                # line = str(line, encoding='ascii')
                ob.append(line.split())
        comments = self.ig.get_comments()
        name = u'REM Name: ' + u'AlOCCF334'
        ob.insert(0, name.split())  # Name is always at first place
        self.assertEqual(comments, ob)

    def testrun_get_firstlast(self):
        files = self.ig.get_gradefiles('./test-data/PFA.gradeserver_all.tgz')
        atoms = self.ig.get_pdbatoms(files[0])
        fl = dbfile.get_first_last_atom(atoms)
        self.assertTupleEqual(fl, ('AL1', 'F36'))

    def testrun_deleted_pdb_file(self):
        with self.assertRaises(SystemExit):
            gdb = dbfile.ParseDB('../userdb.txt')
            ImportGRADE('./PFA.gradeserver_all_2.tgz', gdb)

    def testrun_get_restaraints(self):
        self.maxDiff = None
        restr = self.ig.get_restraints()
        # to generate the test file:
        # with open('test.txt', 'wb+') as file:
        #    for line in restr:
        #        file.write(' '.join(line)+'\n')
        filename = './grade_restraints.txt'
        tst = []
        with open(filename) as test_file:
            for line in test_file:
                tst.append(line.split())
        self.assertListEqual(restr, tst)

    def testrun_get_pdbatoms(self):
        # pdblines = []
        with open('./grade-PFA.pdb') as pdb_file:
            pdblines = pdb_file.readlines()
            pdbatoms = self.ig.get_pdbatoms(pdblines)
            self.assertListEqual(['AL1', 'AL', '9.463', '-3.351', '3.397'], pdbatoms[0])


class DSRParseTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.res_file = 'dsrparse.res'
        testresfile = './dsrparse.res'
        self.rl = ResList(testresfile)
        self.reslist = self.rl.get_res_list()
        self.dsrp = DSRParser(self.reslist)

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
        self.dsrp = DSRParser(self.reslist)

    def testrun_find_dsr_command(self):
        self.maxDiff = None
        line = self.dsrp.find_dsr_command(line=True)
        # self.assertEqual(num, 264)
        string = 'rem dsr put oc(cf3)3 with o1 c1 c2 c3 c4 on O1_3 c1_3 q6 Q4 q7 resi cf3  PART 2 occ -31 dfix\n'
        self.assertEqual(string, line)


class ExportTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.invert = False
        self.gdb = dbfile.ParseDB('../dsr_db.txt')
        self.export_clip = 'benzene'
        self.resgood = ['TITL toluene\n', 'REM This file was exported by DSR version {}\n'.format(VERSION),
                        'REM Name: Toluene, C7H8\nREM Source: CCDC CESLUJ\n',
                        'CELL 0.71073    11.246   14.123   27.184   90.000  100.079   90.000\n',
                        'ZERR    1.00   0.000    0.000    0.000    0.000    0.000    0.000\n', 'LATT  -1\n', 'SFAC C\n',
                        'UNIT 1 \n',
                        'REM  RESIDUE: TOL\n',
                        'REM Sum formula: C7 \n',
                        'WGHT  0.1\n',
                        'FVAR  1\n',
                        'rem Restraints from DSR database:\n',
                        'SADI C2 C3 C3 C4 C4 C5 C5 C6 C6 C7 C7 C2\nSADI 0.04 C2 C6 C2 C4 C7 C5 C3 C7 C4 C6 C3 C5\nDFIX 1.51 C1 C2\nSADI 0.04 C1 C7 C1 C3\nFLAT C1 > C7\nSIMU C1 > C7\nRIGU C1 > C7\n',
                        'rem Restraints from atom connectivities:\n',
                        ['DFIX 1.3922 C3   C2  \n',
                         'DFIX 1.3775 C3   C4  \n',
                         'DFIX 1.5058 C2   C1  \n',
                         'DFIX 1.3946 C2   C7  \n',
                         'DFIX 1.3802 C7   C6  \n',
                         'DFIX 1.3814 C6   C5  \n',
                         'DFIX 1.3897 C5   C4  \n',
                         'DANG 2.5246 C1   C3  \n',
                         'DANG 2.5243 C1   C7  \n',
                         'DANG 2.4183 C2   C4  \n',
                         'DANG 2.4124 C2   C6  \n',
                         'DANG 2.3878 C3   C7  \n',
                         'DANG 2.3909 C3   C5  \n',
                         'DANG 2.3961 C4   C6  \n',
                         'DANG 2.3967 C5   C7  \n',
                         'FLAT C2 C7 C6 C1\n',
                         'FLAT C6 C5 C4 C3\n'],
                        'rem end of restraints\n',
                        '\n',
                        ['C1   1     0.34810   0.50619   0.44851   11.0   0.04\n',
                         'C2   1     0.37174   0.58816   0.41613   11.0   0.04\n',
                         'C3   1     0.27706   0.63878   0.38821   11.0   0.04\n',
                         'C4   1     0.29758   0.71355   0.35825   11.0   0.04\n',
                         'C5   1     0.41548   0.73951   0.35559   11.0   0.04\n',
                         'C6   1     0.51068   0.69033   0.38312   11.0   0.04\n',
                         'C7   1     0.48938   0.61536   0.41297   11.0   0.04\n'],
                        '\nHKLF 0\nEND\n']
        self.resgoodall = ['TITL toluene\n', 'REM This file was exported by DSR version {}\n'.format(VERSION),
                           'REM Name: Toluene, C7H8\nREM Source: CCDC CESLUJ\n',
                           'CELL 0.71073    11.246   14.123   27.184   90.000  100.079   90.000\n',
                           'ZERR    1.00   0.000    0.000    0.000    0.000    0.000    0.000\n',
                           'LATT  -1\n', 'SFAC C\n',
                           'UNIT 1 \n',
                           'REM  RESIDUE: TOL\n',
                           'REM Sum formula: C7 \n',
                           'WGHT  0.1\n',
                           'FVAR  1\n',
                           '\n',
                           ['C1   1     0.34810   0.50619   0.44851   11.0   0.04\n',
                            'C2   1     0.37174   0.58816   0.41613   11.0   0.04\n',
                            'C3   1     0.27706   0.63878   0.38821   11.0   0.04\n',
                            'C4   1     0.29758   0.71355   0.35825   11.0   0.04\n',
                            'C5   1     0.41548   0.73951   0.35559   11.0   0.04\n',
                            'C6   1     0.51068   0.69033   0.38312   11.0   0.04\n',
                            'C7   1     0.48938   0.61536   0.41297   11.0   0.04\n'],
                           '\nHKLF 0\nEND\n']

    def testrun_format_calced_coords(self):
        export = Export(self.gdb, invert=False)
        bigcell = export.format_calced_coords([1, 1, 1, 90, 90, 90], "benzene")[0]
        smallcell = export.format_calced_coords([2, 1, 1, 90, 90, 90], "benzene")[0]
        cell1 = [50, 50, 50, 90, 90, 90]
        cell2 = [2, 1, 1, 90, 90, 90]
        self.assertListEqual(bigcell, cell1)
        self.assertListEqual(smallcell, cell2)

    def testrun_export_to_clip(self):
        """
        Exports the current fragment to the clipboard.
        """
        gdb = dbfile.ParseDB('../dsr_db.txt')
        export = Export(gdb)
        #        with self.assertRaises(SystemExit):
        self.assertTrue(export.export_to_clip('benzene'))


class ResListEditTest(unittest.TestCase):
    def setUp(self):
        self.res_file = './p21c.res'
        self.rl = ResList(self.res_file)
        self.res_list = self.rl.get_res_list()
        self.fa = FindAtoms(self.res_list)
        self.rle = ResListEdit(self.res_list, self.fa)

    def testrun_find_fvarlines(self):
        lines = self.rle.find_fvarlines()
        self.assertListEqual([17], lines)
        self.assertNotEqual(None, lines)
        self.assertNotEqual(False, lines)


class ResidueTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.res_file = './p21c.res'
        self.rl = ResList(self.res_file)
        self.res_list = self.rl.get_res_list()
        self.find_atoms = FindAtoms(self.res_list)
        self.rle = ResListEdit(self.res_list, self.find_atoms)
        self.dsrp = DSRParser(self.res_list)
        fragment = 'ch2cl2'
        invert = False
        self.head2 = ['DFIX 1.771 0.01 CL1 C1 CL2 C1',
                      'DFIX 2.916 0.03 CL1 CL2',
                      'SIMU CL1 > C1',
                      'RIGU CL1 > C1']
        self.gdb = dbfile.ParseDB('../dsr_db.txt')
        self.residue_class = self.gdb.get_resi(fragment)
        self.fragline = self.gdb.get_startline(fragment)  # full string of FRAG line
        self.dbatoms = self.gdb.get_atoms(fragment)  # only the atoms of the dbentry as list
        self.dbhead = self.gdb.get_startline(fragment)
        self.resi = Resi(self.dsrp, self.residue_class, self.find_atoms)

    def testrun_get_resinumber(self):
        self.assertEqual(self.resi.get_resinumber, '4')
        self.assertNotEqual(self.resi.get_resinumber, 'False')
        self.assertEqual(self.resi.get_residue_class, 'CF3')
        self.assertNotEqual(self.resi.get_residue_class, 'CCF3')

    def testrun_format_restraints(self):
        resihead1 = self.resi.format_restraints(self.head2)
        resihead2 = ['DFIX_CF3 1.771 0.01 CL1 C1 CL2 C1',
                     'DFIX_CF3 2.916 0.03 CL1 CL2',
                     'SIMU_CF3 CL1 > C1',
                     'RIGU_CF3 CL1 > C1']
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
        self.assertEqual(self.resi.get_resi_syntax(['CF3', '2', '3']), class_number_alias)
        self.assertEqual(self.resi.get_resi_syntax(['CF3']), only_class)

    # @unittest.skip('foo')
    def testrun_build_up_residue(self):
        dsr_dict1 = {'target': ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7'], 'fragment': 'OC(CF3)3',
                     'occupancy': '-31', 'source': ['O1', 'C1', 'C2', 'C3', 'C4'],
                     'resi': ['CF13', '2'], 'command': 'PUT', 'dfix': False, 'part': '2'}
        dsr_dict2 = {'target': ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7'], 'fragment': 'OC(CF3)3',
                     'occupancy': '-31', 'source': ['O1', 'C1', 'C2', 'C3', 'C4'],
                     'resi': ['CF23', '5'], 'command': 'PUT', 'dfix': False, 'part': '2'}
        self.dsrp.dsr_dict = dsr_dict1
        resi1 = Resi(self.dsrp, self.residue_class, self.find_atoms)
        self.dsrp.dsr_dict = dsr_dict2
        resi2 = Resi(self.dsrp, self.residue_class, self.find_atoms)
        residue1 = {'alias': None, 'class': 'CF13', 'number': '4'}
        residue2 = {'alias': None, 'class': 'CF23', 'number': '5'}
        self.assertDictEqual(resi1.build_up_residue(), residue1)
        self.assertDictEqual(resi2.build_up_residue(), residue2)

    # @unittest.skip('foo')
    def testrun_make_resihead(self):
        testhead = ['RESI TST', 'SADI 0.02 C1 C2 C1 C3 C1 C4',
                    'SADI 0.02 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4',
                    'SADI 0.04 C2 C3 C3 C4 C2 C4', 'SIMU O1 > F9', 'RIGU O1 > F9']
        dsr_dict1 = {'target': ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7'],
                     'fragment': 'OC(CF3)3',
                     'occupancy': '-31',
                     'source': ['O1', 'C1', 'C2', 'C3', 'C4'],
                     'resi': ['CF13'],  # only class given
                     'command': 'PUT',
                     'dfix': False,
                     'part': '2'}
        dsr_dict2 = {'target': ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7'],
                     'fragment': 'OC(CF3)3',
                     'occupancy': '-31',
                     'source': ['O1', 'C1', 'C2', 'C3', 'C4'],
                     'resi': False,  # no resi given
                     'command': 'PUT',
                     'dfix': False,
                     'part': '2'}
        dsr_dict3 = {'target': ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7'],
                     'fragment': 'OC(CF3)3',
                     'occupancy': '-31',
                     'source': ['O1', 'C1', 'C2', 'C3', 'C4'],
                     'resi': '',  # only resi active, class from db
                     'command': 'PUT',
                     'dfix': False,
                     'part': '2'}
        dsr_dict4 = {'target': ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7'],
                     'fragment': 'OC(CF3)3',
                     'occupancy': '-31',
                     'source': ['O1', 'C1', 'C2', 'C3', 'C4'],
                     'resi': ['tes1', '8'],  # class and number given
                     'command': 'PUT',
                     'dfix': False,
                     'part': '2'}
        dsr_dict5 = {'target': ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7'],
                     'fragment': 'OC(CF3)3',
                     'occupancy': '-31',
                     'source': ['O1', 'C1', 'C2', 'C3', 'C4'],
                     'resi': ['89'],  # class and number given
                     'command': 'PUT',
                     'dfix': False,
                     'part': '2'}
        dsr_dict6 = {'target': ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7'],
                     'fragment': 'OC(CF3)3',
                     'occupancy': '-31',
                     'source': ['O1', 'C1', 'C2', 'C3', 'C4'],
                     'resi': ['tes1', '8', '10'],  # class, number and alias given
                     'command': 'PUT',
                     'dfix': False,
                     'part': '2'}
        self.dsrp.dsr_dict = dsr_dict1
        resi1 = Resi(self.dsrp, 'CF13', self.find_atoms)
        self.dsrp.dsr_dict = dsr_dict2
        resi2 = Resi(self.dsrp, 'TEST', self.find_atoms)
        self.dsrp.dsr_dict = dsr_dict3
        resi3 = Resi(self.dsrp, 'TEST', self.find_atoms)
        self.dsrp.dsr_dict = dsr_dict4
        resi4 = Resi(self.dsrp, 'TES1', self.find_atoms)
        self.dsrp.dsr_dict = dsr_dict5
        resi5 = Resi(self.dsrp, 'TES1', self.find_atoms)
        self.dsrp.dsr_dict = dsr_dict6
        resi6 = Resi(self.dsrp, 'TES1', self.find_atoms)


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


if __name__ == "__main__":
    unittest.main()
