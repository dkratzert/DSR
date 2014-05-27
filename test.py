#/usr/bin/env python
#-*- encoding: utf-8 -*-


import sys, os
import unittest
from atomhandling import get_atomtypes, FindAtoms
from resfile import ResList



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
        self.assertFalse(get_atomtypes(self.bad_dbatoms1))
    def testrun_neg2(self):
        self.assertFalse(get_atomtypes(self.bad_dbatoms2))

class reslistTest(unittest.TestCase):
    def setUp(self):
        #print ("setUp reslistTest executed!")
        self.atom1 = 'F10   4   -0.362398   0.278516   0.447770  11.00000   0.02302   0.04023 =\n'
        self.atom2 = '0.02897   0.00131  -0.01216   0.00374\n'
        self.res_file = 'p21c.res'
        self.res_list = ResList(self.res_file)
        self.reslist =  self.res_list.get_res_list()
        self.fa = FindAtoms(self.reslist)
        self.resi_str = 'RESI TOL 3'
        self.resi_list = ['RESI', 'TOL', '3']
        self.resi_list_class = ['RESI', 'TOL']
        self.resi_list_num = ['RESI', '3']
    
    def testrun_pos(self):
        self.assertEqual(self.fa.get_atom(self.atom1), ['F10', '4', '-0.362398', '0.278516', '0.447770'])
        self.assertIsNone(self.fa.get_atom(self.atom2))
    
    def testrun_resinum(self):
        #print(self.fa.get_resinum(self.resi_str))
        self.assertEqual(self.fa.get_resinum(self.resi_str)['number'], '3')
        self.assertEqual(self.fa.get_resinum(self.resi_list)['number'], '3')
        self.assertIsNone(self.fa.get_resinum(self.resi_list_class)['number'])
        self.assertIs(self.fa.get_resinum(self.resi_list_num)['number'], '3')
    
    def testrun_resinum_class(self):
        self.assertEqual(self.fa.get_resinum(self.resi_str)['class'], 'TOL')
        self.assertEqual(self.fa.get_resinum(self.resi_list)['class'], 'TOL')
        self.assertIs(self.fa.get_resinum(self.resi_list_class)['class'], 'TOL')
        self.assertIsNone(self.fa.get_resinum(self.resi_list_num)['class'])

#get_atoms_resiclass
#collect_residues

    def testrun_get_atoms_resinumber(self):
        self.assertEqual(self.fa.get_atoms_resinumber('C1_1'), '1')
        self.assertEqual(self.fa.get_atoms_resinumber('C1_1b'), '1')
        self.assertEqual(self.fa.get_atoms_resinumber('C1_b'), '0')
        self.assertEqual(self.fa.get_atoms_resinumber('C1_'), '0')
        self.assertIsNot(self.fa.get_atoms_resinumber('C1_1'), 1)
        self.assertIsNot(self.fa.get_atoms_resinumber('C1'), 0)
        self.assertEqual(self.fa.get_atoms_resinumber('C1_0'), '0')
        self.assertEqual(self.fa.get_atoms_resinumber('C1'), '0')



if __name__ == "__main__": 
    unittest.main()
