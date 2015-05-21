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

DFIX 1.328 C22 F1A C22 F2A C22 F3A  C22 F4A C22 F5A C22 F6A
DFIX 2.125 F1A F5A F5A F3A F3A F4A F4A F2A F2A F6A F6A F1A
SADI 0.1 C19 F1A C19 F2A C19 F3A  C19 F4A C19 F5A C19 F6A

DFIX 1.328 C1 F1B C1 F2B C1 F3B  C1 F4B C1 F5B C1 F6B
DFIX 2.125 F1B F5B F5B F3B F3B F4B F4B F2B F2B F6B F6B F1B
SADI 0.1 C2 F1B C2 F2B C2 F3B  C2 F4B C2 F5B C2 F6B

basefilename
reslist

'''
import resfile
from dbfile import global_DB
from atomhandling import FindAtoms, SfacTable, Elem_2_Sfac, NumberScheme
from dsrparse import DSR_Parser
from options import OptionsParser
from refine import ShelxlRefine
from restraints import ListFile
from elements import ELEMENTS
from misc import atomic_distance, matrix_mult
from math import sin, cos
import sys
from resfile import ResList


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
        self.make_afix(afixnum='130', linenumber=atomlinenumber[0])

    def make_afix(self, afixnum, linenumber):
        '''
        create an afix to build a CF3 or CH3 group
        :param afixnum: afix number
        :type afixnum: string
        TODO: 
        - find next unused free variable to refine parts
        - refine
        '''
        num_130 = NumberScheme(self.reslist, ['F1', 'F2', 'F3'], False)
        num_120 = NumberScheme(self.reslist, ['F1', 'F2', 'F3', 'F4', 'F5', 'F6'], False)
        # returns also the atom names if residue is active
        numberscheme_120 = num_120.get_fragment_number_scheme()
        numberscheme_130 = num_130.get_fragment_number_scheme()
        sfac = self.e2s.elem_2_sfac('F')
        afix_130 = ['\nAFIX {0}',
                numberscheme_130[0]+' {1} 0 0 0 11 0.04',
                numberscheme_130[1]+' {1} 0 0 0 11 0.04',
                numberscheme_130[2]+' {1} 0 0 0 11 0.04',
                'AFIX 0\n']
        afix_120 = ['\nAFIX {0}',
                'PART 1',
                numberscheme_120[0]+' {1} 0 0 0 11 0.04',
                numberscheme_120[1]+' {1} 0 0 0 11 0.04',
                numberscheme_120[2]+' {1} 0 0 0 11 0.04',
                'PART 2',
                numberscheme_120[3]+' {1} 0 0 0 11 0.04',
                numberscheme_120[4]+' {1} 0 0 0 11 0.04',
                numberscheme_120[5]+' {1} 0 0 0 11 0.04',
                'PART 0',
                'AFIX 0\n']
        afix_130 = '\n'.join(afix_130)
        afix_120 = '\n'.join(afix_120)
        if str(afixnum) == '130':
            afix = afix_130
        elif str(afixnum) == '120':
            afix = afix_120
        else:
            print('Only CF3 groups implemented yet.')
        atomline = self.reslist[linenumber].split()
        if atomline[-1] == '=':
            self.reslist[linenumber] = '{:5.4s}{:4.2s}{:>10.8s} {:>10.8s} {:>10.8s}  {:8.6s}  0.04'.format(*atomline)
            self.reslist[linenumber+1] = '' 
        self.reslist[linenumber] = self.reslist[linenumber]+afix.format(afixnum, sfac) 
        #print(self.reslist[linenumber])

    
    def rotate_fluorine_atom(self, coords):
        '''
        rotates the coordinates of a fluorine atom in the cone 
        around the terminal atom
        '''
        theta = 115
        phi = 10
        Rot = [[cos(phi), sin(phi), 0], 
               [-cos(theta)*sin(phi), cos(theta)*cos(phi), sin(theta)], 
               [sin(theta)*sin(phi), -sin(theta)*cos(phi), cos(theta)]]
        rotated = matrix_mult(Rot, coords)
        return rotated

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

    def make_refine_cycle(rl, reslist):
        shx = ShelxlRefine(reslist, basefilename, find_atoms)
        acta_lines = shx.remove_acta_card()
        shx.set_refinement_cycles('0')
        rl.write_resfile(reslist, '.ins')
        #sys.exit()
        shx.run_shelxl()
        lf = ListFile(basefilename)
        lst_file = lf.read_lst_file()
        shx.check_refinement_results(lst_file)
        rl = ResList(res_file)
        reslist = rl.get_res_list()
        shx.restore_acta_card(acta_lines)
        shx.set_refinement_cycles('8')
        rl.write_resfile(reslist, '.res')

    ####################################################
    cf3 = CF3(rle, find_atoms, reslist, fragment, sfac_table)
    cf3.cf3('C22')
    make_refine_cycle(rl, reslist)

   
    #dsrp.find_dsr_command()


    print('finished...')
    ####################################################
    
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
        