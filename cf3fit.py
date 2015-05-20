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
from atomhandling import FindAtoms, SfacTable, Elem_2_Sfac
from dsrparse import DSR_Parser
from options import OptionsParser
from refine import ShelxlRefine
from restraints import ListFile
from misc import atomic_distance
from elements import ELEMENTS
import sys


class CF3(object):
    '''
    a class to create cf3 groups at terminal atoms
    '''

    def __init__(self, rle, fa, reslist, fragment, sfac_table):
        '''
        Constructor
        '''
        self.fa = fa
        self.fragment = fragment
        self.reslist = reslist
        self.e2s = Elem_2_Sfac(sfac_table)
        atoms = fa.atoms_as_residues
        atomlist = []
        # atomlist: ['C1', ['x', 'y', 'z'], linenumber, class, part, element, sfac_number, residue_num]
        for i in atoms:
            for y in atoms[i]:
                if y[0][0] == 'Q':
                    continue
                atomlist.append(y+[i])
        self.atomlist = atomlist
        self.cell = rle.get_cell()
    
    
    def find_bonded_fluorine(self, atom, extra_param=0.16, element='F'):
        '''
        find fluorine atoms that are boneded to atom
        returns ['C1', ['x', 'y', 'z'], linenumber, class, part, element, sfac_number, residue_num]
        '''
        found_atoms = []
        atcoord = self.fa.get_atomcoordinates([atom])
        cr = ELEMENTS['C'].covrad
        fr = ELEMENTS[element].covrad
        for i in self.atomlist:
            if not i[5] == element:
                continue
            d = atomic_distance(i[1], atcoord[atom], self.cell)
            if d <= (fr+cr)+extra_param and d > (cr or fr):
                found_atoms.append(i)
        return found_atoms
    
    
    def delete_bound_fluorine(self, bound_atoms):
        '''
        deletes fluorine atoms bound to atom
        :param bound_atoms:
        :type bound_atoms:
        '''
        atoms = [str(i[0])+'_'+str(i[7]) for i in bound_atoms]
        target_lines = self.fa.get_atom_line_numbers(atoms)
        for i in target_lines:
            i = int(i)
            rle.remove_line(i, rem=False, remove=False, frontspace=True)
    
    def cf3(self, atom):
        '''
        create CF3 group on atom 
        '''
        atomlinenumber = self.fa.get_atom_line_numbers([atom])
        print(atomlinenumber, '#####', atom)
        found = self.find_bonded_fluorine(atom)
        for i in found:
            print(i[0]+'_'+i[7])
        self.delete_bound_fluorine(found)
        self.make_afix(afixnum='137', linenumber=atomlinenumber[0])

    def make_afix(self, afixnum, linenumber):
        '''
        create an afix to build a CF3 or CH3 group
        :param afixnum: afix number
        :type afixnum: string
        '''
        sfac = self.e2s.elem_2_sfac('F')
        afix_137 = ['\nAFIX {0}',
                'Fxx1 {1} 0 0 0 11 0.04',
                'Fxx2 {1} 0 0 0 11 0.04',
                'Fxx3 {1} 0 0 0 11 0.04',
                'AFIX 0\n']
        afix_137 = '\n'.join(afix_137)
        if str(afixnum) == '137':
            afix = afix_137
        else:
            print('Only CF3 groups implemented yet.')
        atomline = self.reslist[linenumber].split()
        if atomline[-1] == '=':
            self.reslist[linenumber] = '{:5.4s}{:4.2s}{:>10.8s} {:>10.8s} {:>10.8s}  {:8.6s}  0.04'.format(*atomline)
            self.reslist[linenumber+1] = '' 
        self.reslist[linenumber] = self.reslist[linenumber]+afix.format(afixnum, sfac)
        #print(self.reslist[linenumber])

if __name__ == '__main__':
    options = OptionsParser()
    #res_file = options.res_file
    #res_file = '/tmp/mlcp57.res' 
    res_file = 'p21n_cf3.res'
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
    ####################################################
    cf3 = CF3(rle, find_atoms, reslist, fragment, sfac_table)
    cf3.cf3('C22')
    ####################################################
    shx = ShelxlRefine(reslist, basefilename, find_atoms)
    acta_lines = shx.remove_acta_card()
    shx.set_refinement_cycles('0')
    rl.write_resfile(reslist, '.ins')
    sys.exit()
    shx.run_shelxl()
    lf = ListFile(basefilename)
    lst_file = lf.read_lst_file()
    shx.check_refinement_results(lst_file)
    shx.restore_acta_card(acta_lines)
    shx.set_refinement_cycles('8')
    
    #cr = ELEMENTS['C'].covrad
    #fr = ELEMENTS['F'].covrad
    #print(cr+fr+0.16)
    
    
    #cf3.cf6('atom')
    # apply to more than one atom:
    #for at in atomlist:
    #    cf3.cf3('atom')
    
    
            #self.conntable = restr.get_conntable_from_atoms(coords, 
        #                               [i[5] for i in atomlist], # types C, N, O
        #                               [i[0]+'_'+i[7] for i in atomlist], # names C1, N2, O1_3
        #                               extra_param=0.16)
        