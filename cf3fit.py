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
from atomhandling import FindAtoms
from dsrparse import DSR_Parser
from options import OptionsParser


class CF3(object):
    '''
    a class to create cf3 groups at terminal atoms
    '''

    def __init__(self, reslist):
        '''
        Constructor
        '''
        print('hello')
        #print(reslist)
        
    
    


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
    cf3 = CF3(reslist)
    
    #cf3.cf3('atom')
    #cf3.cf6('atom')
    # apply to more than one atom:
    #for at in atomlist:
    #    cf3.cf3('atom')