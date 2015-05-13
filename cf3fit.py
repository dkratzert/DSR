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


class CF3(object):
    '''
    a class to create cf3 groups at terminal atoms
    '''

    def __init__(self, basefilename):
        '''
        Constructor
        '''
        
        
    
    


if '__name__' == '__main__':
    cf3 = CF3()
    
    cf3.cf3('atom')
    cf3.cf6('atom')
    for at in atomlist:
        cf3.cf3('atom')