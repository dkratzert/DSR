import unittest

from dsr_shelx import dbfile
from dsr_shelx.atomhandling import FindAtoms, check_source_target
from dsr_shelx.dsrparse import DSRParser
from dsr_shelx.resfile import ResList, ResListEdit
from dsr_shelx.resi import Resi


class ResidueTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.res_file = './tests/p21c.res'
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
        self.gdb = dbfile.ParseDB('./src/dsr_shelx/dsr_db.txt')
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
        empty_dict = {'alias': '', 'class': '', 'number': ''}
        only_number = {'alias': '', 'class': '', 'number': '2'}
        class_number = {'alias': '', 'class': 'CF3', 'number': '2'}
        only_class = {'alias': '', 'class': 'CF3', 'number': ''}
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


class collect_residuesTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.res_file = './tests/collect_resi.res'
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
