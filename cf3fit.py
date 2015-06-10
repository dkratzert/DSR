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

- wrap restraints after 79 characters

'''
import resfile
from dbfile import global_DB
from atomhandling import FindAtoms, SfacTable, Elem_2_Sfac, NumberScheme
from dsrparse import DSR_Parser
from options import OptionsParser
from refine import ShelxlRefine
from restraints import ListFile
from elements import ELEMENTS
from misc import atomic_distance, frac_to_cart, cart_to_frac,\
    id_generator, shift, remove_partsymbol, find_multi_lines, wrap_headlines
from math import sin, cos, radians, sqrt
import sys
from resfile import ResList, ResListEdit
import mpmath as mp 


dfixr_130 = ['DFIX 1.328 Z F1 Z F2 Z F3 \n', 
             'DFIX 2.125 F1 F2 F2 F3 F3 F1 \n',
             'SADI 0.1 Y F1 Y F2 Y F3 \n',
             'RIGU Y Z F1 F2 F3 ']
dfixr_120 = ['DFIX 1.328 Z F1 Z F2 Z F3  Z F4 Z F5 Z F6 \n', 
             'DFIX 2.125 F1 F2 F2 F3 F3 F1  F4 F5 F5 F6 F6 F4 \n',
             'SADI 0.1 Y F1 Y F2 Y F3  Z F4 Z F5 Z F6 \n',
             'RIGU Y Z F1 F2 F3 \n',
             'EADP F1 F4 \n',
             'EADP F6 F3 \n',
             'EADP F2 F5 ']
dfixr_cf9 = ['SUMP 1 0.0001 1 {0} 1 {1} 1 {2}', 
             'DFIX 1.328 Z F1 Z F2 Z F3  Z F4 Z F5 Z F6 \n', 
             'DFIX 2.125 F1 F2 F2 F3 F3 F1  F4 F5 F5 F6 F6 F4 \n',
             'SADI 0.1 Y F1 Y F2 Y F3  Z F4 Z F5 Z F6 \n',
             'RIGU Y Z F1 F2 F3 \n',
             'EADP F1 F4 \n',
             'EADP F6 F3 \n',
             'EADP F2 F5 ']

sadir_130 = ['SADI 0.02 Z F1 Z F2 Z F3 \n',
             'SADI 0.04 F1 F2 F2 F3 F3 F1 \n',
             'SADI 0.1 Z F1 Z F2 Z F3 \n',
             'RIGU Y Z F1 F2 F3 ']
sadir_120 = ['SADI 0.02 Z F1 Z F2 Z F3  Z F4 Z F5 Z F6 \n',
             'SADI 0.04 F1 F2 F2 F3 F3 F1  F4 F5 F5 F6 F6 F4 \n',
             'SADI 0.1 Z F1 Z F2 Z F3  Z F4 Z F5 Z F6 \n',
             'RIGU Y Z F1 F2 F3 F4 F5 F6 \n',
             'EADP F1 F4 \n',
             'EADP F6 F3 \n',
             'EADP F2 F5 ']
sadir_cf9 = ['SUMP 1 0.0001 1 {0} 1 {1} 1 {2}',
             'SADI 0.02 Z F1 Z F2 Z F3  Z F4 Z F5 Z F6  Z F7 Z F8 Z F9\n',
             'SADI 0.04 F1 F2 F2 F3 F3 F1  F4 F5 F5 F6 F6 F4  F7 F8 F8 F9 F9 F7 \n',
             'SADI 0.1 Z F1 Z F2 Z F3  Z F4 Z F5 Z F6  Z F7 Z F8 Z F9 \n',
             'RIGU Y Z F1 > F9 \n']


class CF3(object):
    '''
    a class to create cf3 groups at terminal atoms
    '''

    def __init__(self, rle, fa, reslist, fragment, sfac_table, basefilename, dsr_dict, resi):
        '''
        Constructor
        '''
        self.resi = resi
        self.rand_id = id_generator(size=7)
        self.fa = fa
        self.rle = rle
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
        self.startm = '\nREM CF3 group made by DSR:\n'
        self.endm = 'REM End of CF3 group made by DSR\n\n'
    
    
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
        for i in bound_atoms:
            i = int(i[2])
            rle.remove_line(i, rem=False, remove=False, frontspace=True)
    
    def cf3(self, afix=130):
        '''
        create CF3 group on atom.
        Either define atom at startup or let dsrparser get the atom name.
        Y-Z-F1/F2/F3
        AFIX 130: CF3 group
        AFIX 120: CF6 group

        :param atom: central atom of the group
        :type atom: string
        :param afix: afix number for the CF3 group
        :type afix: string
        '''
        afix=str(afix)
        if afix == '130':
            print('Generating CF3-Group at {}.'.format(self.dsr_dict['target'][0]))
            restr = sadir_130
            if self.dsr_dict['dfix']:
                restr = dfixr_130
        if afix == '120':
            print('Generating twofold disordered CF3-Group at {}.'.format(self.dsr_dict['target'][0]))
            restr = sadir_120
            if self.dsr_dict['dfix']:
                restr = dfixr_120
        if self.resi.get_residue_class:
            restr = self.resi.format_restraints(restr)
            restr = [i+'\n' for i in restr]
        atom = self.dsr_dict['target'][0]
        if len(self.dsr_dict['target']) > 1:
            print('Using only first target atom {}.'.format(self.dsr_dict['target'][0]))
        atomline = self.fa.get_atom_line_numbers([atom])
        found = self.find_bonded_fluorine(atom)
        for i in found:
            print('Deleting ' + i[0] + '_' + i[7] + ' from '+atom)
        self.delete_bound_fluorine(found)
        fatoms = self.make_afix(afixnum=afix, linenumber=atomline[0])
        self.do_refine_cycle(self.rl, self.reslist)
        # this is the bond around the CF3 group rotates
        Y, Z = self.lf.get_bondvector(atom)
        Y = remove_partsymbol(Y)
        Z = remove_partsymbol(Z)
        if afix == '130':
            F1, F2, F3 = fatoms
            replacelist = (('Z', Z), ('Y', Y), ('F1', F1), ('F2', F2), ('F3', F3))
        if afix == '120':
            F1, F2, F3, F4, F5, F6 = fatoms
            replacelist = (('Z', Z), ('Y', Y), ('F1', F1), ('F2', F2), ('F3', F3),
                            ('F4', F4), ('F5', F5), ('F6', F6))
        for old, new in (replacelist):
            restr = [i.replace(old, new) for i in restr]
        restr = wrap_headlines(restr, 77)
        self.reslist = self.rl.get_res_list()
        # get position for the fluorine atoms:
        atomline = self.fa.get_atom_line_numbers([atom])
        # add restraints to reslist:
        self.reslist[atomline[0]] = self.reslist[atomline[0]]+self.startm+''.join(restr)
        regex = r'.*{}'.format(self.rand_id)
        id_lines = find_multi_lines(self.reslist, regex)
        for line in id_lines:
            self.reslist[line] = ' '.join(self.reslist[line].split()[1:-1])+'\n'
        self.rl.write_resfile(self.reslist, '.res')
        return fatoms

    def cf9(self):
        '''
        create disorderd CF3 group on three positions.
        either define atom at startup or let dsrparser get the atom name.
        
        :param atom: central atom of the group
        :type atom: string
        '''
        print('Generating threefold disordered CF3-Group at {}.'.format(self.dsr_dict['target'][0]))
        atom = self.dsr_dict['target'][0]
        if self.dsr_dict['dfix']:
            restr = dfixr_cf9
        else:
            restr = sadir_cf9 
        if self.resi.get_residue_class:
            restr = self.resi.format_restraints(restr)
            restr = [i+'\n' for i in restr]
        if len(self.dsr_dict['target']) > 1:
            print('Using only first target atom {}.'.format(self.dsr_dict['target'][0]))
        atomline = self.fa.get_atom_line_numbers([atom])
        found = self.find_bonded_fluorine(atom)
        for i in found:
            print('Deleting ' + i[0] + '_' + i[7] + ' from '+atom)
        self.delete_bound_fluorine(found)
        reslist_copy = self.reslist[:]
        fatoms = self.make_afix(afixnum=130, linenumber=atomline[0])
        self.do_refine_cycle(self.rl, self.reslist)
        #self.reslist = self.rl.get_res_list()
        # this is the bond around the CF3 group rotates
        Y, Z = self.lf.get_bondvector(atom)
        Y = remove_partsymbol(Y)
        Z = remove_partsymbol(Z)
        numberedatoms = ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9']
        nums = NumberScheme(self.reslist, numberedatoms, False)
        F1, F2, F3, F4, F5, F6, F7, F8, F9  = nums.get_fragment_number_scheme()
        #start_f_coord = self.fa.get_atomcoordinates([fatoms[0]]).values()[0]
        start_f_coord = self.lf.get_single_coordinate(fatoms[0])
        replacelist = (('Z', Z), ('Y', Y), ('F1', F1), ('F2', F2), ('F3', F3),
                            ('F4', F4), ('F5', F5), ('F6', F6), ('F7', F7), ('F8', F8), ('F9', F9))
        for old, new in (replacelist):
            restr = [i.replace(old, new) for i in restr]
        # add restraints to reslist:
        self.reslist = reslist_copy
        if self.dsr_dict['occupancy']:
            occ = self.dsr_dict['occupancy']
        else:
            occ = str((self.rle.get_fvar_count())*10+1+30)
        fcount = self.rle.get_fvar_count()
        fvar = self.rle.set_free_variables(occ) 
        atomline = self.fa.get_atom_line_numbers([atom])[0]
        atoms_cf9 = ['PART 1 {1}1', 
                    F1+' {0}   {4}   11.0  0.04',
                    F2+' {0}   {5}   11.0  0.04',
                    F3+' {0}   {6}   11.0  0.04',
                    'PART 2 {2}1',
                    F4+' {0}   {7}   11.0  0.04',
                    F5+' {0}   {8}   11.0  0.04',
                    F6+' {0}   {9}   11.0  0.04',
                    'PART 3 {3}1',
                    F7+' {0}   {10}   11.0  0.04',
                    F8+' {0}   {11}   11.0  0.04',
                    F9+' {0}   {12}   11.0  0.04',                    
                    'PART 0\n\n']
        restr = wrap_headlines(restr, 77)
        restr = ''.join(restr).format(fcount+1, fcount+2, fcount+3)
        self.reslist[atomline] = self.reslist[atomline]+self.startm+restr
        at1 = self.lf.get_single_coordinate(Y)
        at2 = self.lf.get_single_coordinate(Z)
        coords = []
        for delta in range(0, 360, 36):
            coord = self.rotate_atom_around_bond(start_f_coord, at1, at2, delta)
            coords.append('{}  {}  {}'.format(*coord))        
        self.reslist[atomline] = self.reslist[atomline]\
                                        + '\n'.join(atoms_cf9).format(self.e2s.elem_2_sfac('F'),
                                                           fcount+1, fcount+2, fcount+3,
                                                           *coords)   
        self.reslist[self.rle.find_fvarlines()[0]] = fvar
        self.rl.write_resfile(self.reslist, '.res')
        return fatoms


    def make_afix(self, afixnum, linenumber):
        '''
        create an afix to build a CF3 or CH3 group
        :param afixnum: afix number
        :type afixnum: string
        '''
        resistr = ''
        resi0 = ''
        occ = 11
        if self.dsr_dict['occupancy']:
            occ = self.dsr_dict['occupancy']
        else:
            occ = str((self.rle.get_fvar_count()+1)*10+1)
        num_130 = NumberScheme(self.reslist, ['F1', 'F2', 'F3'], False)
        num_120 = NumberScheme(self.reslist, ['F1', 'F2', 'F3', 'F4', 'F5', 'F6'], False)
        # returns also the atom names if residue is active
        numberscheme_120 = num_120.get_fragment_number_scheme()
        numberscheme_130 = num_130.get_fragment_number_scheme()
        if self.resi.get_residue_class:
            resiclass = self.resi.get_residue_class
            resinum = self.resi.get_resinumber
            resistr = 'RESI '+resiclass+' '+resinum
            resi0 = 'RESI 0'
            numberscheme_130 = ['F1', 'F2', 'F3']
            numberscheme_120 = ['F1', 'F2', 'F3', 'F4', 'F5', 'F6']
        sfac = self.e2s.elem_2_sfac('F')
        # CF3:
        afix_130 = ['\n'+resistr,
                    'AFIX {0}', # AFIX 120 or 130
                    #'REM AFIX made by DSR: {3}',    
                    numberscheme_130[0]+' {1} 0 0 0 11.0  0.04',
                    numberscheme_130[1]+' {1} 0 0 0 11.0  0.04',
                    numberscheme_130[2]+' {1} 0 0 0 11.0  0.04',
                    #'REM end of AFIX by DSR {3}', # insert ID and later change it to the PART usw.
                    'AFIX 0',
                    resi0+'\n']
        # CF6:
        afix_120 = ['\n'+resistr,
                    'AFIX {0}',
                    '\nREM PART 1 !{3}', 
                    numberscheme_120[0]+' {1} 0 0 0   {2}  0.04',
                    numberscheme_120[1]+' {1} 0 0 0   {2}  0.04',
                    numberscheme_120[2]+' {1} 0 0 0   {2}  0.04',
                    'REM PART 2 !{3}',
                    numberscheme_120[3]+' {1} 0 0 0  -{2}  0.04',
                    numberscheme_120[4]+' {1} 0 0 0  -{2}  0.04',
                    numberscheme_120[5]+' {1} 0 0 0  -{2}  0.04',
                    'REM PART 0 !{3}',
                    'AFIX 0',
                    resi0,
                    self.endm]
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
        self.reslist[linenumber] = self.reslist[linenumber]\
                                        + afix.format(afixnum, sfac, occ, 
                                            #this ID is to recognize this line later:
                                            self.rand_id ) 
        if str(afixnum) == '120':
            return numberscheme_120
        if str(afixnum) == '130':
            return numberscheme_130

        

    
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
        #self.do_refine_cycle(self.rl, self.reslist)
        print(fluorine_names[0])
        ratom = self.lf.get_single_coordinate(fluorine_names[0])
        print('ratom:', ratom)
        #ratom = self.get_coordinates_of_first_atom(fluorine_names)
        self.lf.read_lst_file()
        bondvec = self.lf.get_bondvector()
        print('bondvec:', bondvec)
        at1 = self.lf.get_single_coordinate(bondvec[0].upper())
        at2 = self.lf.get_single_coordinate(bondvec[1].upper())
        names = ['F{}'.format(i) for i in range(1, 25)]
        num_thor = NumberScheme(self.reslist, names, False)
        names = num_thor.get_fragment_number_scheme()
        for delta in range(0, 360, 15):
            coord = self.rotate_atom_around_bond(ratom, at1, at2, delta)
            coords.append(coord)
        #print('PART 0')
        diffden = self.lf.get_difference_density(averaged=False)
        print(diffden)
        rotation = self.lf.get_degree_of_highest_peak()
        print(rotation)
        n = (rotation/15)+4 #+1
        diffden = shift(diffden, n)
        print(len(diffden))
        print('\n')
        print
        print('AFIX 6')
        # make restraints:
        sad12 = []
        sad13 = []
        ff = []
        flast = False
        for i, num, co, dif in zip(names, [i for i in range(1, 25)], coords, diffden):
            #print('PART {}'.format(num))
            if not flast:
                flast = 'F24'
            ff.append((flast, i))
            flast = i
            sad12.append(('C22', i))
            sad13.append(('C19', i))
            print('{}  3  {:0<8.6}  {:0<8.6}  {:0<8.6}  {:<8.6}  -1.2\
            '.format(i, co[0], co[1], co[2], 10.0+(dif/2760.0)))
            # with free variables:
            #print('{}  3  {:0<8.6}  {:0<8.6}  {:0<8.6}  {:<8.6}  261\
            #'.format(i, co[0], co[1], co[2], 10.0*num+1+10))
        from misc import flatten
        print('SADI '+' '.join(flatten(sad12)))
        print('SADI '+' '.join(flatten(sad13)))
        print('SADI '+' '.join(flatten(ff)))
        d = 0
        for dif in diffden:
            d += dif/2760.0
        #print('PART 0')
        print('AFIX 0')
        print
        print(d)

        
        

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
        v = [round(i, 6) for i in v]
        return v 
    
if __name__ == '__main__':
    from resi import Resi
    options = OptionsParser()
    #res_file = options.res_file
    #res_file = '/tmp/mlcp57.res' 
    res_file = 'p21n_cf3.res'
    invert = options.invert
    basefilename = resfile.filename_wo_ending(res_file)
    basefilename = 'p21n_cf3'
    gdb = global_DB(invert)
    rl = resfile.ResList(res_file)
    reslist = rl.get_res_list()
    find_atoms = FindAtoms(reslist)
    rle = ResListEdit(reslist, find_atoms)
    fvarlines = rle.find_fvarlines()
    dsrp = DSR_Parser(reslist, rle)
    dsr_dict = dsrp.get_dsr_dict
    dsr_line_number = dsrp.find_dsr_command(line=False)
    fvarlines = rle.find_fvarlines()
    if dsrp.occupancy:
        rle.set_free_variables(dsrp.occupancy)
    fragment = dsrp.fragment.lower()
    sf = SfacTable(reslist, ['C', 'F', 'F', 'F'])
    sfac_table = sf.set_sfac_table() 

    resi = Resi(reslist, dsr_dict, dbhead='RESI CF3', db_residue_string='CF3', find_atoms=find_atoms)


    ####################################################
    
    
    cf3 = CF3(rle, find_atoms, reslist, fragment, sfac_table, basefilename, dsr_dict, resi)
    
    if fragment == 'cf3':
        cf3.cf9()
    if fragment == 'cf6':
        cf3.cf3('120')
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
        