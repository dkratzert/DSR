import unittest
from pathlib import Path

import pytest

from dsr_shelx import dbfile
from dsr_shelx.export import Export
from dsr_shelx.version import VERSION

TESTS_DIR = Path(__file__).parent


class ExportTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.invert = False
        self.gdb = dbfile.ParseDB('./src/dsr_shelx/dsr_db.txt')
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
        bigcell = export.expand_calced_cell([1, 1, 1, 90, 90, 90], self.gdb.get_atoms(self.export_clip))[0]
        smallcell = export.expand_calced_cell([2, 1, 1, 90, 90, 90], self.gdb.get_atoms(self.export_clip))[0]
        cell1 = [50, 50, 50, 90, 90, 90]
        cell2 = [2, 1, 1, 90, 90, 90]
        self.assertListEqual(bigcell, cell1)
        self.assertListEqual(smallcell, cell2)

    @unittest.skip('This can fail because of user rights')
    def testrun_export_to_clip(self):
        """
        Exports the current fragment to the clipboard.
        """
        gdb = dbfile.ParseDB('./src/dsr_shelx/dsr_db.txt')
        export = Export(gdb)
        self.assertTrue(export.export_to_clip('benzene'))


def test_write_res_file(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Export(dbfile.ParseDB(str(TESTS_DIR.parent / 'src/dsr_shelx/dsr_db.txt')), invert=True).write_res_file('Benzene')
    text = (tmp_path / 'benzene.res').read_text(encoding='latin1')
    for item in ('TITL benzene', 'SFAC C', 'REM Sum formula: C6', 'DFIX 1.', 'FLAT ', 'HKLF 0'):
        assert item in text
    assert len([line for line in text.splitlines() if line.startswith('C')]) == 6 + 1  # atoms and CELL


def test_export_all_fragments(tmp_path, monkeypatch):
    entry = ('<{0}>\nREM Name: {0}\nRESI {0}\nSADI C1 C2 C2 C3\nFRAG 17 1 1 1 90 90 90\n'
             'C1 1 0.0 0.0 0.0\nC2 1 1.5 0.0 0.0\nC3 1 3.0 0.0 0.0\n</{0}>\n')
    (tmp_path / 'db.txt').write_text(entry.format('fraga') + entry.format('fragb'))
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit):
        Export(dbfile.ParseDB(str(tmp_path / 'db.txt'))).export_all_fragments()
    assert (tmp_path / 'fraga.res').is_file() and (tmp_path / 'fragb.res').is_file()


def test_export_to_clip(monkeypatch):
    clipboard = []
    monkeypatch.setattr('pyperclip.setcb', clipboard.append)
    assert Export(dbfile.ParseDB(str(TESTS_DIR.parent / 'src/dsr_shelx/dsr_db.txt'))).export_to_clip('benzene')
    assert clipboard[0].startswith('FRAG') and clipboard[0].endswith('FEND')
    assert clipboard[0].count('\nC') == 6
