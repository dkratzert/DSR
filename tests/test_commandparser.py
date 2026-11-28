import unittest

import pytest

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


def _parse(command):
    return DSRParser([command + '\n', 'HKLF 4\n'])


def test_occupancy_without_leading_zero():
    dsrp = _parse('rem dsr put toluene with C1 C2 C3 on Q1 Q2 Q3 part 1 occ .5')
    assert dsrp.occupancy == '.5'


@pytest.mark.parametrize('occ', ['abc', '-2a1', 'nan', 'inf'])
def test_non_numeric_occupancy_exits_with_message(occ, capsys):
    with pytest.raises(SystemExit):
        _parse('rem dsr put toluene with C1 C2 C3 on Q1 Q2 Q3 part 1 occ ' + occ)
    assert 'Occupancy without numerical value' in capsys.readouterr().out


@pytest.mark.parametrize('command', ['rem dsr put', 'rem dsr'])
def test_incomplete_dsr_command_exits_with_message(command, capsys):
    with pytest.raises(SystemExit):
        _parse(command)
    assert 'Incomplete DSR command' in capsys.readouterr().out
