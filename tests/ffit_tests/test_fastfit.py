import os
import sys
import unittest
from pathlib import Path

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from src.dsr.misc import remove_file, copy_file

print(sys.version)


# @unittest.skip("skipping dsr_complete_runs_Test ")
class dsr_complete_runs_ffit_Test(unittest.TestCase):
    def setUpClass(_=None) -> None:
        pass
        """os.system(f'rsync -rumv --delete-after {"/Users/daniel/Documents/GitHub/DSR"} {"/Applications/"}')
        os.system('rm /Applications/DSR/dsr')
        os.system('cp /Users/daniel/Documents/GitHub/DSR/setup/dsr-mac /Applications/DSR/dsr')
        os.system('chmod a+x /Applications/DSR/dsr')"""

    def setUp(self):
        self.maxDiff = None
        # remove this to view the results:
        # self.dsr = '/Applications/DSR/dsr'
        # self.dsr = 'D:\Programme\DSR\dsr'
        self.dsr = f'python3 src/dsr/dsr.py'
        self.prefix='./tests/ffit_tests'
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
        # 13 -r resi cf3 PART 2 occ 501
        # 14 -r resi cf3 PART 2
        # 15 -r only PUT
        # 16 -r RESI
        # 17 -r resi occ -31
        # 18 occ -31
        # 19 occ -31 without FVAR
        # 20 put CF6 on C22 split (AFIX 120)
        # 21 put CF3 on C22 (AFIX 130)
        # 22 put CF9 on C1
        #  This test will fail with fractional source coordinates:
        # 23 -r -target x y z  PART -1 OCC 10.5 RESI TOL
        # 24 correctly restore res file after SHELXL failure
        # 25 REM DSR PUT CH2CL2 WITH CL1 C1 CL2 ON C01P Q1 C01M PART 1 OCC 71 RESI CCL2
        #    Test for problems with upper/lower case source atom names
        # 26 regular -r dsr run with
        #         dfix PART 2 occ -31     dfix and part without resi
        # 27 -g (rigid) -re part 2 occ -31
        # 28 PART -1 OCC 10.5 DFIX -> negative part and dfix

    def dsr_runtest(self, nummer=99, parameter='-r', external_file='', hkl=None,
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
        if not remlines:
            remlines = []
        d = []
        c = []
        # parameter = '-noffit ' + parameter
        print(f'{nummer} ' * 10, 'start:')
        print(f'## Running test in: {Path(".").resolve()}')
        if hkl:
            copy_file(f'{self.prefix}/{hkl}.hkl', f'{self.prefix}/{nummer}a.hkl')
        copy_file(f'{self.prefix}/{nummer}.res', f'{self.prefix}/{nummer}a.res')
        os.system(f'{self.dsr} {parameter} {self.prefix}/{nummer}a.res')
        with open(f'{self.prefix}/{nummer}a.{ending}') as txt:
            a = txt.readlines()[limit_start:limit_end]
            a = [x.strip(' \n\r') for x in a]
        with open(f'{self.prefix}/{nummer}-erg.{ending}') as txt2:
            b = txt2.readlines()[limit_start:limit_end]
            b = [x.strip(' \n\r') for x in b]
        if external_file:
            with open(f'{self.prefix}/{external_file}.dfix') as ext:
                c = ext.readlines()
            with open(f'{self.prefix}/{external_file}-erg.dfix') as ext2:
                d = ext2.readlines()
        for line in remlines:
            a[line] = ''
            b[line] = ''
        print(f'{nummer} test:')
        print("parameter:", parameter)
        if hkl:
            remove_file(f'{self.prefix}/{nummer}a.hkl')
        remove_file(f'{self.prefix}/{nummer}a.fcf')
        remove_file(f'{self.prefix}/{nummer}.fcf')
        remove_file(f'{self.prefix}/{nummer}.2fcf')
        remove_file(f'{self.prefix}/{nummer}a.lst')
        # a = remove_whitespace(a)
        # b = remove_whitespace(b)
        self.assertEqual('\n'.join(b), '\n'.join(a))
        if external_file:
            self.assertEqual(d, c)
        print('{} '.format(nummer) * 10, "ende")
        remove_file(f'{self.prefix}/{nummer}a.ins')
        remove_file(f'{self.prefix}/{nummer}a.res')
        remove_file(f'{self.prefix}/{external_file}.dfix')

    # @unittest.skip(" skipping1 ")
    def testrun_run1(self):
        """
        regular dsr run with
        resi cf3 PART 2 occ -31
        """
        self.maxDiff = None
        self.dsr_runtest(1, '-r', remlines=[])

    # @unittest.skip(" skipping2 ")
    def testrun_run2(self):
        """
        regular -r dsr run with
        resi cf3 dfix =
            PART 2 occ -31
        """
        self.maxDiff = None
        self.dsr_runtest(2, '-r', remlines=[])

    # @unittest.skip(" skipping3 ")
    def testrun_run3(self):
        """
        regular run with:
         occ -31 PART 2
        """
        self.dsr_runtest(3, '-r', remlines=[])

    # @unittest.skip(" skipping4 ")
    def testrun_run4(self):
        """
        external -re restraints with:
        resi cf3 =
            PART 2 occ -31
        """
        self.dsr_runtest(4, '-re', external_file='dsr_CCF3_6_4a', remlines=[])

    # @unittest.skip(" skipping5 ")
    def testrun_run5(self):
        """
        dsr -re  with:
        REM dsr put oc(CF3)3 with o1 c1 c2 C3 on O1_4 C1_4 Q6 Q7 resi Ccf3 part 2 occ =
            -31 dfix
        """
        self.dsr_runtest(5, '-re', external_file='dsr_CCF3_6_5a_dfx', remlines=[])

    # @unittest.skip(" skipping6 ")
    def testrun_run6(self):
        """
        -re   PART 2 occ -31
        """
        self.dsr_runtest(6, '-re', external_file='dsr_1_6a', remlines=[])

    # @unittest.skip(" skipping 7")
    def testrun_run7(self):
        """
        rigid

        """
        self.dsr_runtest(7, '-g -r', remlines=[])

    # @unittest.skip(" skipping 8")
    def testrun_8(self):
        """
        dsr -s tol

        """
        print('8 ' * 10, 'start')
        os.system(f"{self.dsr} -s tol > search.txt")
        with open('search.txt') as txt:
            se = txt.readlines()
        good = False
        print(''.join(se))
        for line in se:
            if line.startswith(" toluene           | Toluene, C7H8"):
                good = True
        if good:
            remove_file('./search.txt')
        self.assertTrue(good, "Search text differs")
        print('8 ' * 10, 'ende')

    # @unittest.skip(" skipping 9")
    def testrun_9(self):
        """
        dsr -e tol

        """
        self.maxDiff = None
        print('9 ' * 10, 'start')
        current_dir = Path(".").resolve()
        print(f'## Running test in: {current_dir}')
        print(f"{self.dsr} -e tolUene")
        os.system(f"{self.dsr} -e tolUene")
        ex = Path('toluene.res').read_text().splitlines(keepends=False)
        ex_erg = Path(f'{self.prefix}/toluene-erg.res').read_text().splitlines(keepends=False)
        del ex[1]  # line with the version number
        del ex_erg[1]
        Path('toluene.res').unlink(missing_ok=True)
        self.assertEqual('\n'.join(ex), '\n'.join(ex_erg))
        print('9 ' * 10, 'ende')

    # @unittest.skip(" skipping 10")
    def testrun_run10(self):
        """
        invert fragment
        -r -t
        """
        # Line 74 has one digit difference in Windows and Mac:
        self.dsr_runtest(10, '-t -r', remlines=[])

    # @unittest.skip(" skipping 11")
    def testrun_run11(self):
        """
        regular dsr run with
        replace resi PART 0
        """
        self.dsr_runtest(11, '-r', remlines=[])

    # @unittest.skip(" skipping 12")
    def testrun_run12(self):
        """
        regular dsr run without fit with
        resi cf3 PART 2 occ -31
        """
        self.dsr_runtest(12, ending='res', parameter='-n -r', remlines=[])

    # @unittest.skip(" skipping 13")
    def testrun_run13(self):
        """
        rem dsr put oc(cf3)3 with o1 c1 c2 c3 c4 on O1_3 c1_3 q6 Q4 q7 resi cf3 =
            PART 2 occ 501
        """
        self.dsr_runtest(13, '-r', remlines=[])

    # @unittest.skip(" skipping 14")
    def testrun_run14(self):
        """
        dsr -r p21c.res
        REM  dsr put oc(CF3)3 with o1 c1 c2 C3 on O1_4 C1_4 Q6 Q7 resi cf3 =
            PART 2
        """
        self.dsr_runtest(14, '-r', remlines=[])

    # @unittest.skip(" skipping 15")
    def testrun_run15(self):
        """
        dsr -r foo
        REM  dsr put oc(CF3)3 with o1 c1 c2 C3 on O1_4 C1_4 Q6 Q7
        """
        self.dsr_runtest(15, '-r', remlines=[])

    # @unittest.skip(" skipping 16")
    def testrun_run16(self):
        """
        dsr -r
        REM  dsr put oc(CF3)3 with o1 c1 c2 C3 on O1_4 C1_4 Q6 Q7 RESI
        """
        self.dsr_runtest(16, '-r', remlines=[])

    # @unittest.skip(" skipping 17")
    def testrun_run17(self):
        """
        dsr -r with:
        REM  dsr put oc(CF3)3 with o1 c1 c2 C3 on O1_4 C1_4 Q6 Q7 resi occ -31
        """
        self.dsr_runtest(17, '-r', remlines=[])

    # @unittest.skip(" skipping 18")
    def testrun_run18(self):
        """
        dsr -r foo
        REM  dsr put oc(CF3)3 with o1 c1 c2 C3 on O1_4 C1_4 Q6 Q7 occ -31
        """
        self.dsr_runtest(18, '-r', remlines=[])

    # @unittest.skip(" skipping 19")
    def testrun_run19(self):
        """
        occ -31 without FVAR

        """
        self.dsr_runtest(19, '-re', external_file='dsr_1_19a', remlines=[])

    # @unittest.skip(" skipping 20")
    def testrun_run20(self):
        """
        rem dsr put CF6 on C22 split
        (AFIX 120)

        """
        self.dsr_runtest(20, '-r', hkl=20, remlines=[149, 150])

    # @unittest.skip(" skipping 21")
    def testrun_run21(self):
        """
        rem dsr put CF3 on C22
        (AFIX 130)
        """
        self.dsr_runtest(21, '-r', hkl=20, remlines=[27, 137, 138])

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
        self.dsr_runtest(23, '-target 1.0005 0.5447 0.5342 0.9314 0.5395 0.5126 0.9995 0.4553 0.4658 -r',
                         hkl=23, remlines=[])

    # @unittest.skip(" skipping 24")
    def testrun_run24(self):
        """
        rem dsr put toluene on C1 C2 C3

        """
        self.dsr_runtest(24, '-r')

    # @unittest.skip(" skipping 25")
    def testrun_run25(self):
        """
        REM DSR PUT CH2CL2 WITH CL1 C1 CL2 ON C01P Q1 C01M PART 1 OCC 71 RESI CCL2
        Test for problems with upper/lower case source atom names
        """
        self.dsr_runtest(25, '-target 0.12127 0.47704 0.70694 0.38298 0.29146 0.52138 0.2512 0.3015 0.6978 -r')

    # @unittest.skip(" skipping26 ")
    def testrun_run26(self):
        """
        regular -r dsr run with
        dfix PART 2 occ -31
        """
        self.dsr_runtest(26, '-r', remlines=[])

    # @unittest.skip(" skipping27 ")
    def testrun_run27(self):
        """
        -re   PART 2 occ -31
        """
        self.dsr_runtest(27, ' -g -re', remlines=[])

    # @unittest.skip(" skipping 28")
    def testrun_run28(self):
        """
        REM DSR PUT TOLUENE WITH C2 C3 C5 ON Q1 Q2 Q1 PART -1 OCC 10.5 DFIX
        1.0005, 0.5447, 0.5342, 0.9314, 0.5395, 0.5126, 0.9995, 0.4553, 0.4658
        """
        self.dsr_runtest(28, '-target 1.0005 0.5447 0.5342 0.9314 0.5395 0.5126 0.9995 0.4553 0.4658 -r', remlines=[])

    # @unittest.skip(" skipping 29")
    def testrun_run29(self):
        """
        REM DSR PUT PENTAFL WITH C1 C5 C3 ON C02Y C04L C03W PART 1 OCC 11 RESI PEFL
        0.99892, 0.65486, 0.64458, 1.03021, 0.55299, 0.63399, 0.89952, 0.57358, 0.69892
        notice: There will be no HAFIX, because the "REM Restraints for Fragment ..." is
                already there. 
        """
        self.dsr_runtest(29, '-target 0.99892 0.65486 0.64458 1.03021 0.55299 0.63399 0.89952 0.57358 0.69892 -r',
                         remlines=[])


def remove_whitespace(mystringlist):
    newlist = []
    for line in mystringlist:
        line = ' '.join(line.split()).strip(' \r\n')
        if line == ' ' or line == '':
            continue
        newlist.append(line)
    return newlist


if __name__ == "__main__":
    unittest.main()
