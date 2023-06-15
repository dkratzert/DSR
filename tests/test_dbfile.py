import unittest

from dsr_shelx import dbfile
from dsr_shelx.atomhandling import rename_restraints_atoms


class TestrenameDBHeadatoms(unittest.TestCase):
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


class dbfileTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.rdb = dbfile.ParseDB('./tests/db1.txt', "./tests/db2.txt")
        self.klein = ['\n', '<DMX>\n', 'REM test\n', 'RESI 3 TST1\n', 'SIMU C1\n', 'FRAG 17 1 1 1 90 90 90\n',
                      'O1  1  -1.3542148   -0.4780990   -0.5279749\n', '</DMX>']

    def testrun_dublicate_tags(self):
        with self.assertRaises(SystemExit):
            rdb = dbfile.ParseDB("./tests/db1_dublicate.TXT")


class globalDB(unittest.TestCase):
    def setUp(self):
        self.maxDiff = 9999
        self.klein = ['RESI ClBE ', 'SADI C1 Cl1 C2 Cl2 ',
                      'SADI 0.04 C6 Cl1 C2 Cl1 C1 Cl2 C3 Cl2 ',
                      'SAME C2 > C6 C1 ', 'FLAT C1 > Cl2 ',
                      'SIMU C1 > Cl2 ', 'RIGU C1 > Cl2']
        self.result = {'dmx': {'startline': 2,
                               'name': 'dmx',
                               'comments': ['REM test'],
                               'atoms': [['O1', 1, -1.3542148, -0.478099, -0.5279749]],
                               'cell': [1.0, 1.0, 1.0, 90.0, 90.0, 90.0],
                               'source': '',
                               'resi': 'TST1',
                               'restraints': ['SIMU C1'],
                               'endline': 8,
                               'hfix': [],
                               'dbname': 'dsr_user_db'},
                       'dme': {'startline': 1,
                               'name': 'dme',
                               'comments': ['REM test'],
                               'atoms': [['O1', 1, -1.3542148, -0.478099, -0.5279749]],
                               'cell': [1.0, 1.0, 1.0, 90.0, 90.0, 90.0],
                               'source': '',
                               'resi': '',
                               'restraints': [],
                               'endline': 5,
                               'hfix': [],
                               'dbname': 'dsr_db'}
                       }

    def testrun_build_db_dict(self):
        gdb = dbfile.ParseDB("tests/db1_klein.TXT", "tests/db2_klein.TXT")
        db = gdb.databases
        self.assertEqual(db['dmx']['startline'], 2)
        self.assertEqual(db['dmx']['name'], 'dmx')  # no name given, so name is fragment tag
        self.assertEqual(db['dmx']['dbname'], 'dsr_user_db')
        self.assertEqual(db['dmx']['restraints'], ['SIMU C1'])
        self.assertEqual(db['dmx']['cell'], [1, 1, 1, 90, 90, 90])
        self.assertEqual(db['dmx']['resi'], 'TST1')
        self.assertEqual(db, self.result)

    def testrun_get_residue_from_head2(self):
        # raises System exit, because residue in db_resinum.TXT is badly defined.
        main_dbpath = "tests/db_resinum.TXT"
        user_dbpath = "tests/db1.TXT"
        with self.assertRaises(SystemExit):
            gdb = dbfile.ParseDB(main_dbpath, user_dbpath)

    def testrun_get_fragment_atoms(self):
        x = -1.154
        z = 0.526
        o1 = ['O1', 1, -1.154, -0.748, 0.526]
        gdb = dbfile.ParseDB("tests/db1.TXT", "tests/db2.TXT")
        atom = gdb.get_atoms('dme')[0]
        self.assertListEqual(o1, atom)
        self.assertEqual(x, atom[2])
        self.assertEqual(z, atom[4])
        self.assertEqual(1, atom[1])
        self.assertEqual('O1', atom[0])

    def testrun_get_fragment_atoms_shortline(self):
        gdb = dbfile.ParseDB("tests/db1_shortline.TXT")
        # db = gdb.build_db_dict()
        atom = gdb.get_atoms('dme-free')
        self.assertEqual(len(atom), 5)
        self.assertEqual(atom[0][0], 'O2')

    def testrun_get_fragment_atoms_inv(self):
        x = -1.154
        z = 0.526
        o1 = ['O1', 1, -1.154, -0.748, 0.526]
        gdb = dbfile.ParseDB("tests/db1.TXT", "tests/db2.TXT")
        atom = gdb.get_atoms('dme')[0]
        self.assertListEqual(o1, atom)
        self.assertEqual(x, atom[2])
        self.assertEqual(z, atom[4])
        self.assertEqual(1, atom[1])
        self.assertEqual('O1', atom[0])

    def testrun_get_fragment_atoms_noatoms(self):
        with self.assertRaises(SystemExit):
            gdb = dbfile.ParseDB("tests/db1_noatoms.TXT", "tests/db2.TXT")
            gdb.get_atoms('dme-free')

    def testrun_get_fragment_atoms_noend(self):
        with self.assertRaises(SystemExit):
            gdb = dbfile.ParseDB("tests/db1_noend.TXT", "tests/db2.TXT")
            gdb.get_atoms('dme-free')

    def testrun_header_consistency(self):
        self.maxDiff = None
        main_dbpath = "./tests/db1_head_inconsistent.TXT"
        user_dbpath = "./tests/db2_klein.TXT"
        with self.assertRaises(SystemExit):
            gdb = dbfile.ParseDB(main_dbpath, user_dbpath)
            fragment = 'dmel'
            gdb.check_db_restraints_consistency(fragment)

    def testrun_header_consistency2(self):
        self.maxDiff = None
        maindb = "./tests/db1_head_inconsistent2.TXT"
        userdb = "./tests/db2_klein.TXT"
        with self.assertRaises(SystemExit):
            gdb = dbfile.ParseDB(maindb, userdb)
            fragment = 'dmem'
            head = gdb[fragment]['restraints']
            atoms = gdb.get_atoms(fragment, invert=True)
            gdb.check_db_restraints_consistency(fragment)

    def testrun_get_resi_from_fragment(self):
        self.maxDiff = None
        gdb = dbfile.ParseDB("tests/comment.TXT", 'tests/db1.txt')
        fragment = 'com1'
        resi = gdb.get_resi(fragment)
        line = gdb.get_startline(fragment)
        self.assertEqual(resi, 'DME')
        self.assertEqual(line, 2)
