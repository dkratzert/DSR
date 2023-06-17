import unittest
from pathlib import Path

from dsr_shelx import dbfile
from dsr_shelx.dbfile import ImportGRADE
from misc import read_cif
from tests import grade2_import_test_result


class ImportGRADE_Test(unittest.TestCase):
    def setUp(self):
        self.gdb = dbfile.ParseDB(maindb_path=Path('src/dsr_shelx/dsr_db.txt').resolve().__str__())
        self.ig = ImportGRADE('./tests/test-data/PFA.gradeserver_all.tgz', self.gdb)
        self.igi = ImportGRADE('./tests/test-data/PFA.gradeserver_all.tgz', self.gdb, invert=True)

    def testrun_get_comments(self):
        self.maxDiff = None
        filename = './tests/grade-comments.dfix'
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
        files = self.ig.get_gradefiles(Path('./tests/test-data/PFA.gradeserver_all.tgz'))
        atoms = self.ig.get_pdbatoms(files.pdb_file)
        fl = dbfile.get_first_last_atom(atoms)
        self.assertTupleEqual(fl, ('AL1', 'F36'))

    def test_get_grade_v2_atoms(self):
        files = self.ig.get_gradefiles(Path('./tests/test-data/pcb.gradeserver_v2_all.tgz'))
        atoms = self.ig.get_pdbatoms(files.pdb_file)
        self.assertListEqual(atoms,
                             [['CL1', 'CL', '3.866', '0.027', '-0.274'],
                              ['C2', 'C', '2.147', '-0.006', '-0.118'],
                              ['C3', 'C', '1.448', '1.168', '0.147'],
                              ['CL4', 'CL', '2.299', '2.660', '0.321'],
                              ['C5', 'C', '0.059', '1.139', '0.273'],
                              ['CL6', 'CL', '-0.810', '2.592', '0.602'],
                              ['C7', 'C', '-0.641', '-0.058', '0.135'],
                              ['C8', 'C', '-2.149', '-0.087', '0.272'],
                              ['O9', 'O', '-2.628', '-0.283', '1.409'],
                              ['O10', 'O1-', '-2.825', '0.089', '-0.775'],
                              ['C11', 'C', '1.456', '-1.205', '-0.258'],
                              ['CL12', 'CL', '2.318', '-2.664', '-0.587'],
                              ['C13', 'C', '0.067', '-1.229', '-0.131'],
                              ['CL14', 'CL', '-0.792', '-2.715', '-0.303']])

    def test_get_grade_v2_atom2(self):
        ig = ImportGRADE('./tests/test-data/pcb.gradeserver_v2_all.tgz', self.gdb)
        imported_entry = ig.bild_grade_db_entry()
        dbentry = ig.compile_dbentry(imported_entry, 'pcb'.upper())
        self.assertEqual(grade2_import_test_result.result_v2, imported_entry)
        self.assertEqual(grade2_import_test_result.result_v2_as_text, dbentry)

    def test_get_grade_v2_cif(self):
        ig = ImportGRADE('./tests/test-data/pcb.gradeserver_v2_all.tgz', self.gdb)
        cif = read_cif(ig.gradefiles.cif_file.splitlines(keepends=False))
        self.assertEqual('foobar', cif.get('_chem_comp.name')[0])

    def testrun_deleted_pdb_file(self):
        with self.assertRaises(SystemExit):
            gdb = dbfile.ParseDB('tests/test-data/userdb.txt')
            ImportGRADE('./tests/PFA.gradeserver_all_2.tgz', gdb)

    def testrun_get_restaraints(self):
        self.maxDiff = None
        restr = self.ig.get_restraints()
        # to generate the test file:
        # with open('test.txt', 'wb+') as file:
        #    for line in restr:
        #        file.write(' '.join(line)+'\n')
        filename = './tests/grade_restraints.txt'
        tst = []
        with open(filename) as test_file:
            for line in test_file:
                tst.append(line.split())
        self.assertListEqual(restr, tst)

    def testrun_get_pdbatoms(self):
        # pdblines = []
        with open('./tests/grade-PFA.pdb') as pdb_file:
            pdblines = pdb_file.read()
            pdbatoms = self.ig.get_pdbatoms(pdblines)
            self.assertListEqual(['AL1', 'AL', '9.463', '-3.351', '3.397'], pdbatoms[0])
