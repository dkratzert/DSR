#-*- encoding: utf-8 -*-
'''
Created on 13.05.2015

@author: Daniel Kratzert
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <daniel.kratzert@ac.uni-freiburg.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#

- force name tags to be longer than three characters 

- remove hydrogen and fluorine atoms from the target atom
- run shelxl with hfix 137 or afix 130 to get the position
- calculate the exact position of the fluorine atoms from the
  values of the shelxl run
- place the atoms with their restraints
- set the part, residue und occupancy
- set the sump and free variable
-

basefilename
reslist

'''
import resfile
from dbfile import global_DB
from atomhandling import FindAtoms, SfacTable
from dsrparse import DSR_Parser
from options import OptionsParser
from refine import ShelxlRefine


class CF3(object):
    '''
    a class to create cf3 groups at terminal atoms
    '''

    def __init__(self, reslist, fa, sfac_table):
        '''
        Constructor
        '''
        print('hello')
        self.fa = fa
        self.reslist = reslist
        self.sfac_table = sfac_table
        
        
    
    def cf3(self, atom):
        '''
        create CF3 group on atom 
        ''' 
        print('removing hydrogens at atom {}\n'.format(atom))
        delcount = self.fa.remove_adjacent_hydrogens(atom, self.sfac_table)
        print('removed:', delcount)
        for i in reslist[103:119]:
            print(i.strip('\n'))
    


if __name__ == '__main__':
    options = OptionsParser()
    #res_file = options.res_file
    res_file = 'p21c.res' 
    invert = options.invert
    basefilename = resfile.filename_wo_ending(res_file)
    gdb = global_DB(invert)
    rl = resfile.ResList(res_file)
    reslist = rl.get_res_list()
    find_atoms = FindAtoms(reslist)
    rle = resfile.ResListEdit(reslist, find_atoms)
    dsrp = DSR_Parser(reslist, rle)
    dsr_dict = dsrp.get_dsr_dict
    fvarlines = rle.find_fvarlines()
    if dsrp.occupancy:
        rle.set_free_variables(dsrp.occupancy, fvarlines)
    fragment = dsrp.fragment.lower()
    sf = SfacTable(reslist, ['C', 'F', 'F', 'F'])
    sfac_table = sf.set_sfac_table() 
    cf3 = CF3(reslist, find_atoms, sfac_table)
    
    shx = ShelxlRefine(reslist, basefilename, find_atoms)
    acta_lines = shx.remove_acta_card()
    shx.set_refinement_cycles('0')
    rl.write_resfile(reslist, '.ins')
    
    cf3.cf3(['C37'])
    
    #cf3.cf6('atom')
    # apply to more than one atom:
    #for at in atomlist:
    #    cf3.cf3('atom')