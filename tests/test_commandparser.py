import unittest

from dsr_shelx.dsrparse import DSRParser
from dsr_shelx.resfile import ResList


class DSRParseTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.res_file = 'tests/dsrparse.res'
        testresfile = 'tests/dsrparse.res'
        self.rl = ResList(testresfile)
        self.reslist = self.rl.get_res_list()
        self.dsrp = DSRParser(self.reslist)

    def testrun_find_dsr_command(self):
        self.maxDiff = None
        num = self.dsrp.dsr_line_number
        self.assertEqual(num, 264)


class DSRParse2Test(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.res_file = 'tests/dsrparse.res'
        testresfile = 'tests/dsrparse.res'
        self.rl = ResList(testresfile)
        self.reslist = self.rl.get_res_list()
        self.dsrp = DSRParser(self.reslist)

    def testrun_find_dsr_command(self):
        self.maxDiff = None
        line = self.dsrp.dsr_command_list
        # self.assertEqual(num, 264)
        string = 'rem dsr put oc(cf3)3 with o1 c1 c2 c3 c4 on O1_3 c1_3 q6 Q4 q7 resi cf3  PART 2 occ -31 dfix\n'
        self.assertEqual(string.upper().split(), line)
