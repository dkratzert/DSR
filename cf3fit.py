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

- force name tags to be longer than three characters?
  Maybe better to reserve CF3, CF6 anf CFx for CF3 groups.  

- remove hydrogen and fluorine atoms from the target atom
- make a copy of hkl and res file
- run shelxl with hfix 137 or afix 130 to get the position
- calculate the exact position of the fluorine atoms from the
  values in the list file. Position of one atom is enough. 
- delete the res/hkl/lst file
- place the atoms with their restraints
- set the part, residue und occupancy
- set the sump and free variable
-

1. generate ordinary cf3 with restraints
2. generate cf6 with ideal staggered conformation
3. generate thorus of isotropic atoms with occupation relative to 
   difference density and inside AFIX 6 or 9

Try this:
--------------------------------------------------
BLOC -1 -1 # -1 for only Uij and occupancy
F1 ... 
F2 ... 
F3 ...
[...]
BLOC 0
with 
SUMP 3 1 1 1 2 1 3 1 4 ...
+ one FVAR for all Uij
------------------------------------------------

DFIX 1.328 C22 F1A C22 F2A C22 F3A  C22 F4A C22 F5A C22 F6A
DFIX 2.125 F1A F5A F5A F3A F3A F4A F4A F2A F2A F6A F6A F1A
SADI 0.1 C19 F1A C19 F2A C19 F3A  C19 F4A C19 F5A C19 F6A

DFIX 1.328 C1 F1B C1 F2B C1 F3B  C1 F4B C1 F5B C1 F6B
DFIX 2.125 F1B F5B F5B F3B F3B F4B F4B F2B F2B F6B F6B F1B
SADI 0.1 C2 F1B C2 F2B C2 F3B  C2 F4B C2 F5B C2 F6B

 Difference electron density (eA^-3x100) at 15 degree intervals for AFIX 130
 group attached to C1_a. The center of the range is eclipsed (cis) to C7_a
 and rotation is clockwise looking down C2_a to C1_a.
   349  237  171  203  358  579  668  504  243   65  -11   41  272  530  563  380  223  146  106  146  316  497  522  451

 After local symmetry averaging:    272   149    88   130   315   536   584   445
'''
import resfile
from dbfile import global_DB
from atomhandling import FindAtoms, SfacTable, Elem_2_Sfac, NumberScheme
from dsrparse import DSR_Parser
from options import OptionsParser
from refine import ShelxlRefine
from restraints import ListFile
from elements import ELEMENTS
from misc import atomic_distance, frac_to_cart, cart_to_frac, copy_file,\
    id_generator, shift, remove_partsymbol
from math import sin, cos, radians, sqrt
import sys
from resfile import ResList
import mpmath as mp
from _collections import deque


class CF3(object):
    '''
    a class to create cf3 groups at terminal atoms
    '''

    def __init__(self, rle, fa, reslist, fragment, sfac_table, basefilename, dsr_dict):
        '''
        Constructor
        '''
        self.rand_id = id_generator(size=7)
        self.fa = fa
        self.dsr_dict = dsr_dict
        self.fragment = fragment
        self.reslist = reslist
        self.rl = resfile.ResList(res_file)
        self.e2s = Elem_2_Sfac(sfac_table)
        atoms = fa.atoms_as_residues
        self.basefilename = basefilename
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
    
    def cf3(self, atom=None):
        '''
        create CF3 group on atom.
        Either define atom at startup or let dsrparser get the atom name.
        Y-Z-F1/F2/F3
        :param atom: central atom of the group
        :type atom: string
        '''
        print('Generating CF3-Group')
        restraints = ['DFIX 1.328 Z F1 Z F2 Z F3 \n', 
                      'DFIX 2.125 F1 F2 F2 F3 F3 F1 \n',
                      'SADI 0.1 Y F1 Y F2 Y F3 \n']
        if not atom:
            atom = self.dsr_dict['target'][0]
        atomlinenumber = self.fa.get_atom_line_numbers([atom])
        #print(atomlinenumber, '#####', atom)
        found = self.find_bonded_fluorine(atom)
        for i in found:
            print('Deleting ' + i[0] + '_' + i[7] + ' from '+atom)
        self.delete_bound_fluorine(found)
        fatoms = self.make_afix(afixnum='130', linenumber=atomlinenumber[0])
        self.do_refine_cycle(self.rl, self.reslist)
        Y, Z = self.lf.get_bondvector()
        Y = remove_partsymbol(Y)
        Z = remove_partsymbol(Z)
        F1, F2, F3 = fatoms
        for old, new in (('Z', Z), ('Y', Y), ('F1', F1), ('F2', F2), ('F3', F3)):
            restraints = [i.replace(old, new) for i in restraints]
        self.reslist = self.rl.get_res_list()
        atomlinenumber = self.fa.get_atom_line_numbers([atom])
        self.reslist[atomlinenumber[0]] = self.reslist[atomlinenumber[0]]+'\n'+''.join(restraints)
        self.rl.write_resfile(self.reslist, '.res')
        return fatoms
        
    def cf6(self, atom=None):
        '''
        create disorderd CF3 group on two positions.
        either define atom at startup or let dsrparser get the atom name.
        
        :param atom: central atom of the group
        :type atom: string
        '''
        print('Generating disordered CF3-Group')
        if not atom:
            atom = self.dsr_dict['target'][0]
        atomlinenumber = self.fa.get_atom_line_numbers([atom])
        #print(atomlinenumber, '#####', atom)
        found = self.find_bonded_fluorine(atom)
        for i in found:
            print('Deleting ' + i[0] + '_' + i[7] + ' from '+atom)
        self.delete_bound_fluorine(found)
        fatoms = self.make_afix(afixnum='120', linenumber=atomlinenumber[0])
        self.do_refine_cycle(rl, reslist)
        return fatoms
    
    def make_cf3_thorus(self, atom=None):
        '''
        Creates a thorus of isotropic fluorine atoms around the central
        atom of a cf3 group. The occupancy is estimated from the 
        residual density values in the lst file.
        
        :param atom: central atom of the cf3 group
        :type atom: string
        '''
        coords = []
        if not atom:
            atom = self.dsr_dict['target'][0]
        # returns the atom names of the fluorine atoms:
        fluorine_names = self.cf3() 
        self.do_refine_cycle(self.rl, self.reslist)
        print(fluorine_names[0]+'_A')
        ratom = self.lf.get_single_coordinate(fluorine_names[0]+'_A')
        print(ratom)
        #ratom = self.get_coordinates_of_first_atom(fluorine_names)
        self.lf.read_lst_file()
        bondvec = self.lf.get_bondvector()
        print(bondvec)
        at1 = self.lf.get_single_coordinate(bondvec[0].upper())
        at2 = self.lf.get_single_coordinate(bondvec[1].upper())
        names = ['F{}'.format(i) for i in range(1, 25)]
        num_thor = NumberScheme(self.reslist, names, False)
        names = num_thor.get_fragment_number_scheme()
        for delta in range(0, 360, 15):
            coord = self.rotate_atom_around_bond(ratom, at1, at2, delta)
            coords.append(coord)
        print('PART 0')
        diffden = self.lf.get_difference_density(averaged=False)
        print(diffden)
        rotation = self.lf.get_degree_of_highest_peak()
        print(rotation)
        n = (rotation/15) #+1
        diffden = shift(diffden, n)
        print(len(diffden))
        print('\n')
        print
        print('AFIX 6')
        for i, num, co, dif in zip(names, [i for i in range(1, 25)], coords, diffden):
            print('PART {}'.format(num))
            print('{}  3  {:0<8.6}  {:0<8.6}  {:0<8.6}  {:<8.6}  -1.2\
            '.format(i, co[0], co[1], co[2], 10.0+(dif/2760.0)))
            # with free variables:
            print('{}  3  {:0<8.6}  {:0<8.6}  {:0<8.6}  {:<8.6}  261\
            '.format(i, co[0], co[1], co[2], 10.0*num+1+10))
        d = 0
        for dif in diffden:
            d += dif/2760.0
        print('PART 0')
        print('AFIX 0')
        print
        print(d)
        #make atoms from coords here    
        
        

    def do_refine_cycle(self, rl, reslist):
        shx = ShelxlRefine(reslist, self.basefilename, find_atoms)
        acta_lines = shx.remove_acta_card()
        shx.set_refinement_cycles('0')
        rl.write_resfile(reslist, '.ins')
        #sys.exit()
        shx.run_shelxl()
        self.lf = ListFile(self.basefilename)
        lst_file = self.lf.read_lst_file()
        shx.check_refinement_results(lst_file)
        rl = ResList(res_file)
        reslist = rl.get_res_list()
        shx.restore_acta_card(acta_lines)
        shx.set_refinement_cycles('8')
        rl.write_resfile(reslist, '.res')
        
        
    def make_afix(self, afixnum, linenumber, atomtype='F'):
        '''
        create an afix to build a CF3 or CH3 group
        :param afixnum: afix number
        :type afixnum: string
        TODO: 
        - find next unused free variable to refine parts
        - for 120, find atoms after refine cycle and add parts
        '''
        occ = 11
        if self.dsr_dict['occupancy']:
            occ = self.dsr_dict['occupancy']
        num_130 = NumberScheme(self.reslist, ['F1', 'F2', 'F3'], False)
        num_120 = NumberScheme(self.reslist, ['F1', 'F2', 'F3', 'F4', 'F5', 'F6'], False)
        # returns also the atom names if residue is active
        numberscheme_120 = num_120.get_fragment_number_scheme()
        numberscheme_130 = num_130.get_fragment_number_scheme()
        sfac = self.e2s.elem_2_sfac('F')
        afix_130 = ['\nAFIX {0}',
                    'REM AFIX made by DSR: {3}',    
                    numberscheme_130[0]+' {1} 0 0 0 {2}  0.04',
                    numberscheme_130[1]+' {1} 0 0 0 {2}  0.04',
                    numberscheme_130[2]+' {1} 0 0 0 {2}  0.04',
                    'REM end of AFIX by DSR {3}', # insert ID and later change it to the PART usw.
                    'AFIX 0\n']
        # TODO: set the occupancy coorectly
        afix_120 = ['\nAFIX {0}',
                    'REM PART 1 {3}', 
                    numberscheme_120[0]+' {1} 0 0 0  {2}  0.04',
                    numberscheme_120[1]+' {1} 0 0 0  {2}  0.04',
                    numberscheme_120[2]+' {1} 0 0 0  {2}  0.04',
                    'REM PART 2 {3}',
                    numberscheme_120[3]+' {1} 0 0 0  {2}  0.04',
                    numberscheme_120[4]+' {1} 0 0 0  {2}  0.04',
                    numberscheme_120[5]+' {1} 0 0 0  {2}  0.04',
                    'REM PART 0 {3}',
                    'AFIX 0\n']
        afix_130 = '\n'.join(afix_130)
        afix_120 = '\n'.join(afix_120)
        if str(afixnum) == '130':
            afix = afix_130
        elif str(afixnum) == '120':
            afix = afix_120
        else:
            print('Only CF3 groups implemented yet.')
            return False
        atomline = self.reslist[linenumber].split()
        # make shure the pivot atom is isotropic:
        if atomline[-1] == '=':
            self.reslist[linenumber] = '{:5.4s}{:4.2s}{:>10.8s} {:>10.8s} {:>10.8s}  {:8.6s}  0.04'.format(*atomline)
            self.reslist[linenumber+1] = '' 
        # insert the afix:
        self.reslist[linenumber] = self.reslist[linenumber]+afix.format(afixnum, 
                                                                        sfac, 
                                                                        occ, 
                                                                        #this ID is to recognize this line later
                                                                        self.rand_id ) 
        if str(afixnum) == '120':
            return numberscheme_120
        if str(afixnum) == '130':
            return numberscheme_130

    
    
    def rotate_atom_around_bond(self, ratom, at1, at2, delta=10):
        '''
        R = T**-1*Rx**-1*Ry**-1*Rz*Ry*Rx*T
        '''
        ratom = frac_to_cart(ratom, self.cell)
        at1 = frac_to_cart(at1, self.cell)
        at2 = frac_to_cart(at2, self.cell)
                
        ratom = mp.matrix(list(ratom)+[1.0])
        delta = radians(delta)
        x0, y0, z0 = at1
        T = mp.matrix(((1, 0, 0, -x0),
                       (0, 1, 0, -y0),
                       (0, 0, 1, -z0),
                       (0, 0, 0,   1)))
        T1 = mp.inverse(T)

        vx = at2[0]-at1[0] 
        vy = at2[1]-at1[1]  # P2 - P1
        vz = at2[2]-at1[2]
        vnorm = sqrt(vx**2+vy**2+vz**2)
        a, b, c = vx/vnorm, vy/vnorm, vz/vnorm
        d = mp.sqrt(b**2+c**2)
        
        sina = b/d # For rotation around alpha
        cosa = c/d #

        Rxa = mp.matrix(((         1,          0,          0,  0),
                         (         0,       cosa,       sina,  0),
                         (         0,      -sina,       cosa,  0),
                         (         0,          0,          0,  1)))
        '''
        Rya = mp.matrix(((      cosa,          0,      -sina,  0),
                         (         0,          1,          0,  0),
                         (      sina,          0,       cosa,  0),
                         (         0,          0,          0,  1)))
        ''' '''
        Rza = mp.matrix(((      cosa,       sina,          0,  0),
                         (     -sina,       cosa,          0,  0),
                         (         0,          0,          1,  0),
                         (         0,          0,          0,  1)))
        '''
        Rxa1 = mp.inverse((Rxa))
        
        cosb = d  #
        sinb = -a # rotation around beta

        Ryb = mp.matrix(((      cosb,          0,      -sinb,  0),
                         (         0,          1,          0,  0),
                         (      sinb,          0,       cosb,  0),
                         (         0,          0,          0,  1)))

        Ryb1 = mp.inverse((Ryb))
                
        sind = mp.sin(delta)
        cosd = mp.cos(delta)

        Rzd = mp.matrix(((      cosd,       sind,          0,  0),
                         (     -sind,       cosd,          0,  0),
                         (         0,          0,          1,  0),
                         (         0,          0,          0,  1)))
        
        R = T1*Rxa*Ryb*Rzd*Ryb1*Rxa1*T
        v = R*ratom
        v = cart_to_frac(v[:3], self.cell)
        return v 
    
if __name__ == '__main__':
    options = OptionsParser()
    #res_file = options.res_file
    #res_file = '/tmp/mlcp57.res' 
    res_file = 'p21n_cf3.res'
    invert = options.invert
    basefilename = resfile.filename_wo_ending(res_file)
    # spaeter das orginaes resfile verenden und nicht mit dsr-tmp:
    copy_file(basefilename+'.res', 'dsr-tmp.res')
    copy_file(basefilename+'.hkl', 'dsr-tmp.hkl')
    res_file = 'dsr-tmp.res'
    basefilename = 'dsr-tmp'
    gdb = global_DB(invert)
    rl = resfile.ResList(res_file)
    reslist = rl.get_res_list()
    find_atoms = FindAtoms(reslist)
    rle = resfile.ResListEdit(reslist, find_atoms)
    dsrp = DSR_Parser(reslist, rle)
    dsr_dict = dsrp.get_dsr_dict
    dsr_line_number = dsrp.find_dsr_command(line=False)
    fvarlines = rle.find_fvarlines()
    if dsrp.occupancy:
        rle.set_free_variables(dsrp.occupancy, fvarlines)
    fragment = dsrp.fragment.lower()
    sf = SfacTable(reslist, ['C', 'F', 'F', 'F'])
    sfac_table = sf.set_sfac_table() 

    


    ####################################################
    
    
    cf3 = CF3(rle, find_atoms, reslist, fragment, sfac_table, basefilename, dsr_dict)
    
    if fragment == 'cf3':
        cf3.cf3()
    if fragment == 'cf6':
        cf3.cf6()
    #cf3.make_cf3_thorus()


   
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
"""

    for i in range(0, 360, 10):
        co = cf3.rotate_atom_around_bond(F1A, C22, C19, delta=i)
        #co = cf3.rotate_fluorine_atom(F1A, C22, C19, alpha=i)
        print('F{}  3  {:0<8.6}  {:0<8.6}  {:0<8.6}  10.08  -1.5'.format(i, *co))
    
    
    
    
    for i in range(0, 360, 10):
        co = cf3.rotate_atom_around_bond(F1A, C22, C19, delta=i)
        #co = cf3.rotate_fluorine_atom(F1A, C22, C19, alpha=i)
        print('F{}  3  {:0<8.6}  {:0<8.6}  {:0<8.6}  10.08  -1.5'.format(i, *co))
"""
        