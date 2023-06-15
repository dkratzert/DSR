import unittest

from dsr_shelx import dbfile, afix
from dsr_shelx.atomhandling import FindAtoms, get_atomtypes, NumberScheme
from dsr_shelx.dsrparse import DSRParser
from dsr_shelx.resfile import ResList


class removeDublicatesAfixTest(unittest.TestCase):
    def setUp(self):
        # self.verbosity = 4
        self.res_file = './tests/collect_resi.res'
        self.res_list = ResList(self.res_file)
        self.reslist = self.res_list.get_res_list()
        self.find_atoms = FindAtoms(self.reslist)
        invert = False
        self.dsrp = DSRParser(self.reslist)
        fragment = 'OC(cf3)3'
        self.gdb = dbfile.ParseDB('./src/dsr_shelx/dsr_db.txt')
        self.dbatoms = self.gdb.get_atoms(fragment)  # only the atoms of the dbentry as list
        self.dbtypes = get_atomtypes(self.dbatoms)
        # self.sf = SfacTable(self.reslist, self.dbtypes)
        # self.sfac_table = self.sf.set_sfac_table()
        self.sfac_table = ['C', 'H', 'N', 'O', 'F']
        self.resi = 'CCF3'  # gdb.get_resi_from_fragment(fragment)
        self.num = NumberScheme(self.reslist, self.dbatoms, self.dsrp)
        self.numberscheme = self.num.get_fragment_number_scheme()
        self.db_testhead = ['SADI_CCF3 C1 C2 C1 C3 C1 C4',
                            'SADI_CCF3 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4 ',
                            'REM test']
        self.db_testhead2 = ['SADI_CF3 C1 C2 C1 C3 C1 C4 ',
                             'SADI_CF3 F1 C2 F2 C2 F3 C2 F4 C3 F5 C3 F6 C3 F7 C4 F8 C4 F9 C4 ']

    def testrun_remove_dublicate_restraints(self):
        newhead = afix.remove_duplicate_restraints(self.reslist, self.db_testhead)
        newhead2 = afix.remove_duplicate_restraints(self.reslist, self.db_testhead2)
        self.assertListEqual(['REM test'], newhead)
        self.assertListEqual(self.db_testhead2, newhead2)
