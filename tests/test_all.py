# /usr/bin/env python
# -*- encoding: utf-8 -*-


import doctest
import unittest

from dsr_shelx import dsr, afix, dsrparse, export, misc, elements, networkx, atomhandling, dbfile, \
    restraints
from dsr_shelx.atomhandling import get_atomtypes, FindAtoms
from dsr_shelx.fit import quatfit
from dsr_shelx.resfile import ResList


class doctestsTest(unittest.TestCase):
    def testrun_doctest(self):
        print('------------- Running Doctests: -------------------')
        for name in [afix, dsrparse, export, misc, elements, atomhandling, networkx.classes.graph, dbfile, restraints,
                     quatfit]:
            failed, attempted = doctest.testmod(name)  # , verbose=True)
            if failed == 0:
                print('passed all {} tests in {}!'.format(attempted, name.__name__))
            else:
                msg = '!!!!!!!!!!!!!!!! {} of {} tests failed in {}  !!!!!!!!!!!!!!!!!!!!!!!!!!!'.format(failed,
                                                                                                         attempted,
                                                                                                         name.__name__)
                self.assertFalse(failed, msg)


def remove_whitespace(mystringlist):
    newlist = []
    for line in mystringlist:
        line = ' '.join(line.split()).strip(' \r\n')
        if line == ' ' or line == '':
            continue
        newlist.append(line)
    return newlist


db_testhead = ['SADI C1 C2 C1 C3 C1 C4',
               'SADI F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4',
               'SADI 0.04 C2 C3 C3 C4 C2 C4', 'SADI 0.04 O1 C2 O1 C3 O1 C4',
               'SADI 0.04 F1 F2 F2 F3 F3 F1 F4 F5 F5 F6 F6 F4 F7 F8 F8 F9 F9 F7',
               'SADI 0.04 F1 C1 F2 C1 F3 C1 F4 C1 F5 C1 F6 C1 F7 C1 F8 C1 F9 C1',
               'SIMU O1 > F9', 'RIGU O1 > F9']

wraphead = ['SADI C1 C2 C1 C3 C1 C4\n',
            'SADI F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4 F4 C3 F5 C3 F6 =\n   C3 F7 C4 F8 C4 F9 C4\n',
            'SADI 0.04 C2 C3 C3 C4 C2 C4\n', 'SADI 0.04 O1 C2 O1 C3 O1 C4\n',
            'SADI 0.04 F1 F2 F2 F3 F3 F1 F4 F5 F5 F6 F6 F4 F7 F8 F8 F9 F9 F7\n',
            'SADI 0.04 F1 C1 F2 C1 F3 C1 F4 C1 F5 C1 F6 C1 F7 C1 F8 C1 F9 C1\n',
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


@unittest.skip("skipping insertAfixTest")
class insertAfixTest(unittest.TestCase):
    def setUp(self):
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


if __name__ == "__main__":
    unittest.main()
