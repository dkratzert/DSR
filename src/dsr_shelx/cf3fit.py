# -*- encoding: utf-8 -*-
"""
Created on 13.05.2015

@author: Daniel Kratzert
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <dkratzert@gmx.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
"""

import string
import sys

from math import radians, sqrt

import mpmath as mpm
from atomhandling import Elem_2_Sfac, NumberScheme
from elements import get_radius_from_element
from misc import atomic_distance, frac_to_cart, cart_to_frac, \
    id_generator, shift, remove_partsymbol, find_multi_lines, wrap_headlines, ufrac_to_ucart, A
from refine import ShelxlRefine
from resfile import ResList
from restraints import ListFile

# Y-Z-F1/F2/F3

dfixr_130 = ['DFIX 1.328 Z F1 Z F2 Z F3 ',
             'DFIX 2.125 F1 F2 F2 F3 F3 F1 ',
             'SADI 0.1   Y F1 Y F2 Y F3 ',
             'RIGU Y Z F1 F2 F3 ',
             'SIMU Y Z F1 F2 F3 ']

dfixr_120 = ['DFIX 1.328 Z F1 Z F2 Z F3  Z F4 Z F5 Z F6 ',
             'DFIX 2.125 F1 F2 F2 F3 F3 F1  F4 F5 F5 F6 F6 F4 ',
             'SADI 0.1   Y F1 Y F2 Y F3  Y F4 Y F5 Y F6 ',
             'RIGU Y Z F1 > F6',
             'SIMU Y Z F1 > F6']

dfixr_120_split = ['DFIX 1.328 ZA F1 ZA F2 ZA F3  ZB F4 ZB F5 ZB F6 ',
                   'DFIX 2.125 F1 F2 F2 F3 F3 F1  F4 F5 F5 F6 F6 F4 ',
                   'SADI 0.1   Y F1 Y F2 Y F3  Y F4 Y F5 Y F6 ',
                   'SADI Y ZA Y ZB',
                   'EADP ZA ZB',
                   'RIGU Y ZA ZB F1 > F6',
                   'SIMU Y ZA ZB F1 > F6']

dfixr_cf9 = ['SUMP 1 0.0001 1 {0} 1 {1} 1 {2}',
             'DFIX 1.328 Z F1 Z F2 Z F3  Z F4 Z F5 Z F6  Z F7 Z F8 Z F9 ',
             'DFIX 2.125 F1 F2 F2 F3 F3 F1  F4 F5 F5 F6 F6 F4  F7 F8 F8 F9 F9 F7 ',
             'SADI 0.1 Y F1 Y F2 Y F3  Y F4 Y F5 Y F6  Y F7 Y F8 Y F9 ',
             'RIGU Y Z F1 > F9',
             'SIMU Y Z F1 > F9']

sadir_130 = ['SADI 0.02 Z F1 Z F2 Z F3 ',
             'SADI 0.04 F1 F2 F2 F3 F3 F1 ',
             'SADI 0.1  Y F1 Y F2 Y F3 ',
             'RIGU Y Z F1 F2 F3 ',
             'SIMU Y Z F1 F2 F3 ']

sadir_120 = ['SADI 0.02 Z F1 Z F2 Z F3  Z F4 Z F5 Z F6 ',
             'SADI 0.04 F1 F2 F2 F3 F3 F1  F4 F5 F5 F6 F6 F4 ',
             'SADI 0.1  Y F1 Y F2 Y F3  Y F4 Y F5 Y F6 ',
             'RIGU Y Z F1 > F6',
             'SIMU Y Z F1 > F6']

sadir_120_split = ['SADI 0.02 ZA F1 ZA F2 ZA F3  ZB F4 ZB F5 ZB F6 ',
                   'SADI 0.04 F1 F2 F2 F3 F3 F1  F4 F5 F5 F6 F6 F4 ',
                   'SADI 0.1  Y F1 Y F2 Y F3  Y F4 Y F5 Y F6 ',
                   'SADI Y ZA Y ZB',
                   'EADP ZA ZB',
                   'RIGU Y ZA ZB F1 > F6',
                   'SIMU Y ZA ZB F1 > F6']

sadir_cf9 = ['SUMP 1 0.0001 1 {0} 1 {1} 1 {2}',
             'SADI 0.02 Z F1 Z F2 Z F3  Z F4 Z F5 Z F6  Z F7 Z F8 Z F9 ',
             'SADI 0.04 F1 F2 F2 F3 F3 F1  F4 F5 F5 F6 F6 F4  F7 F8 F8 F9 F9 F7 ',
             'SADI 0.1  Y F1 Y F2 Y F3  Y F4 Y F5 Y F6  Y F7 Y F8 Y F9 ',
             'RIGU Y Z F1 > F9',
             'SIMU Y Z F1 > F9']


class CF3(object):
    '''
    a class to create cf3 groups at terminal atoms
    '''

    def __init__(self, rle, fa, reslist, fragment, sfac_table, basefilename, dsrp, resi, res_file, options):
        """
        Constructor
        """
        self.options = options
        self.resi = resi
        self.rand_id = id_generator(size=7)
        self.rand_id_part_neg = id_generator(size=7)
        self.fa = fa
        self.rle = rle
        self.dsrp = dsrp
        self.fragment = fragment
        self.reslist = reslist
        self.res_file = res_file
        self.rl = ResList(res_file)
        self.e2s = Elem_2_Sfac(sfac_table)
        atoms = fa.atoms_as_residues()
        self.basefilename = basefilename
        atomlist = []
        # atomlist: ['C1', ['x', 'y', 'z'], linenumber, class, part, element,
        #                                           sfac_number, residue_num]
        for i in atoms:
            for y in atoms[i]:
                if y[0][0] == 'Q':
                    continue
                atomlist.append(y + [i])
        self.atomlist = atomlist
        self.cell = rle.get_cell()
        self.startm = '\nREM CF3 group made by DSR:\n'
        self.endm = 'REM End of CF3 group made by DSR\n'

    def find_bonded_fluorine(self, atom, extra_param=0.16, element='F'):
        """
        find fluorine atoms that are boneded to atom
        returns ['C1', ['x', 'y', 'z'], linenumber, class, part, element, sfac_number, residue_num]
        """
        found_atoms = []
        atcoord = self.fa.get_atomcoordinates([atom])
        cr = get_radius_from_element('C')
        fr = get_radius_from_element(element)
        for i in self.atomlist:
            if not i[5] == element:
                continue
            d = atomic_distance(i[1], atcoord[atom], self.cell)
            if d <= (fr + cr) + extra_param and d > (cr or fr):
                found_atoms.append(i)
        return found_atoms

    def delete_bound_fluorine(self, bound_atoms):
        """
        deletes fluorine atoms bound to atom
        :param bound_atoms:
        :type bound_atoms:
        """
        for i in bound_atoms:
            i = int(i[2])
            self.rle.remove_line(i, rem=False, remove=True, frontspace=False)

    def format_cf3_restraints(self, afix, restr, atom, fatoms, splitatoms=False):
        """
        replaces the dummy atom names in the restraint lists with the real names
        :param afix: string of the afix number
        :param restr: restraints
        :param atom: pivot atom
        :param fatoms: fluorine atoms
        """
        Y, Z = self.lf.get_bondvector(atom)
        Y = remove_partsymbol(Y)
        Z = remove_partsymbol(Z)
        replacelist = ()
        if afix == '130':
            F1, F2, F3 = fatoms
            replacelist = (('Z', Z), ('Y', Y), ('F1', F1), ('F2', F2), ('F3', F3))
        if afix == '120' and splitatoms:
            F1, F2, F3, F4, F5, F6 = fatoms
            ZA, ZB = splitatoms
            replacelist = (('ZA', ZA), ('ZB', ZB), ('Y', Y), ('F1', F1), ('F2', F2),
                           ('F3', F3), ('F4', F4), ('F5', F5), ('F6', F6))
        if afix == '120' and not splitatoms:
            F1, F2, F3, F4, F5, F6 = fatoms
            replacelist = (('Z', Z), ('Y', Y), ('F1', F1), ('F2', F2), ('F3', F3),
                           ('F4', F4), ('F5', F5), ('F6', F6))
        # replace dummy atoms in restraint list with real atom names:
        for old, new in replacelist:
            restr = [i.replace(old, new) for i in restr]
        restr = wrap_headlines(restr, 77)
        return restr

    def make_pivot_isotropic(self, linenumber: int):
        """
        make sure the pivot atom of a cf3 group is isotropic
        :param linenumber: line number (index) in self.reslist of the pivot atom
        :type linenumber: integer
        :return Uij values, coordinates: U and xyz values of the pivot atom as lists
        """
        atomline = self.reslist[linenumber].split()
        coords = [float(atomline[2]), float(atomline[3]), float(atomline[4])]
        if atomline[-1] == '=':
            nextline = self.reslist[linenumber + 1].split()
            try:
                U11, U22 = atomline[6], atomline[7]
                U33, U23, U13, U12 = nextline[0], nextline[1], nextline[2], nextline[3]
            except Exception:
                # In this case we have a U value missing
                print('*** Incomplete Uij values. Atom split not possible! ***')
                self.dsrp.split = False
                return [[], coords]
            self.reslist[linenumber] = '{:5.4s}{:4.2s}{:>10.8s} {:>10.8s} {:>10.8s}  {:8.6s}  0.04'.format(*atomline)
            self.reslist[linenumber + 1] = ''
            return [[float(U11), float(U22), float(U33), float(U23), float(U13), float(U12)], coords]
        else:
            # atom is already isotropic, nothing to do...
            if self.dsrp.split:
                print('Pivot atom is isotropic. Atom split not possible!')
                self.dsrp.split = False
            return [[], coords]

    def add_chars(self, atom, alphabet):
        """
        add chars to an atom name until the name is unique
        :param atom:
        """
        for char in alphabet:
            if atom + char not in [i[0] for i in self.atomlist]:
                del alphabet[0]
                return atom + char

    def prepare_cf3(self):
        """
        prepares some things before CF3 group is generated
        :return: list
        """
        targetatom = self.dsrp.target[0]
        if targetatom.startswith("Q"):
            print("*** Only carbon atoms allowed for CF3-groups! ***")
            sys.exit()
        if len(self.dsrp.target) > 1:
            print('Using only first target atom {}.'.format(targetatom))
        try:
            atomline = self.fa.get_atom_line_numbers([targetatom])[0]
        except IndexError:
            print("\n*** Atom {} not found ***\n".format(targetatom))
            sys.exit()
        if '_' in targetatom:
            print('\n*** Sorry, can not create a CF3 group inside a residue! '
                  '\nThis would damage the residue. ***')
            sys.exit()
        found = self.find_bonded_fluorine(targetatom)
        for i in found:
            print('Deleting {}_{} from {}'.format(i[0], i[7], targetatom))
        self.delete_bound_fluorine(found)
        return atomline

    def cf3(self, afix='130'):
        """
        create CF3 group on atom.
        Either define atom at startup or let dsrparser get the atom name.
        Y-Z-F1/F2/F3
        AFIX 130: CF3 group
        AFIX 120: CF6 group

        :param afix: afix number for the CF3 group
        :type afix: string
        """
        # the pivot atom of the CF3 group:
        atomline = self.prepare_cf3()
        afix = str(afix)
        splitat1 = ''
        splitat2 = ''
        axes = []
        targetatom = self.dsrp.target[0]
        # c_coords: coordinates of the pivot atom
        uvals, c_coords = self.make_pivot_isotropic(atomline)
        restr = ['']
        if afix == '130':
            print('Generating CF3-Group at {}.'.format(targetatom))
            restr = sadir_130
            if self.dsrp.dfix:
                restr = dfixr_130
        if afix == '120':
            print('Generating twofold disordered CF3-Group at {}.'.format(targetatom))
            restr = sadir_120
            if self.dsrp.split:
                restr = sadir_120_split
            if self.dsrp.dfix:
                restr = dfixr_120
                if self.dsrp.split:
                    restr = dfixr_120_split
        if afix == '120' and self.dsrp.split and uvals:
            num = NumberScheme(self.reslist, [targetatom], self.dsrp)
            if len(targetatom) < 4:
                # in this case it is possible to add a character
                alphabet = [i for i in string.ascii_uppercase]
                splitat1 = self.add_chars(targetatom, alphabet)
                splitat2 = self.add_chars(targetatom, alphabet)
            else:
                splitat1 = num.get_fragment_number_scheme()[0]
                splitat2 = num.get_fragment_number_scheme(extranames=[splitat1])[0]
            splitatoms = [splitat1, splitat2]
            axes = calc_ellipsoid_axes(c_coords, uvals, self.cell)
        else:
            splitatoms = False
        # The fluorine atoms are generated here:
        fatoms = self.make_afix(afixnum=int(afix), linenumber=atomline)
        if not fatoms:
            return False
        self.do_refine_cycle(self.rl, self.reslist)
        # this is essential
        self.reslist = self.rl.get_res_list()
        # this is the bond around the CF3 group rotates
        restr = self.format_cf3_restraints(afix, restr, targetatom, fatoms, splitatoms)
        # get position for the fluorine atoms and make sure the reslist is the newest:
        self.fa._reslist = self.reslist
        atomline = self.fa.get_atom_line_numbers([targetatom])[0]
        # add restraints to reslist:
        restr = ''.join(restr)
        if not self.dsrp.split:
            self.reslist[atomline] = self.reslist[atomline] + self.startm + restr
        else:
            self.reslist[atomline] = ''
            self.reslist[atomline + 1] = self.reslist[atomline] + self.startm + restr
        id_lines = find_multi_lines(self.reslist, r'.*{}'.format(self.rand_id))
        neg_part_id_lines = find_multi_lines(self.reslist, r'.*{}'.format(self.rand_id_part_neg))
        # replace dummy PART with real part definition and C-atom coords with split coords
        if self.dsrp.split and afix == '120':
            at1 = '{:<5s} {:<3} {:>9.6f}   {:>9.6f}   {:>9.6f}   {:>8.4f}     0.04\n' \
                .format(splitat1, self.e2s.elem_2_sfac('C'), axes[0][0],
                        axes[0][1], axes[0][2], float(self.dsrp.occupancy))
            at2 = '{:<5s} {:<3} {:>9.6f}   {:>9.6f}   {:>9.6f}   {:>8.4f}     0.04\n' \
                .format(splitat2, self.e2s.elem_2_sfac('C'), axes[1][0],
                        axes[1][1], axes[1][2], -float(self.dsrp.occupancy))
            dummy = ''
            splb = False
            # put all together to build up the atoms:
            for line, splatom in zip(id_lines, [at1, at2, dummy]):
                # restraints should never be placed in this reslist[line]:
                self.reslist[line] = ' '.join(self.reslist[line].split()[1:3]) + '\n'
                self.reslist[line] = self.reslist[line] + splatom
                # exchange the first three Fluorine atoms with coordinates from one
                # disorder direction and the last three with the other direction:
                for num, fline in enumerate(self.reslist[line + 1:line + 1 + 3]):
                    fline = fline.split()
                    if not splatom:
                        continue
                    [fsplit_a, fsplit_b] = calc_ellipsoid_axes(fline[2:5], uvals, self.cell)
                    resline = self.reslist[line + 1 + num].split()
                    resline[5] = float(resline[5])
                    resline[2:5] = fsplit_a
                    if splb:
                        resline[2:5] = fsplit_b
                    if num == 2 and not splb:
                        # switch to fsplit_b. I can not decide on num, because that goes two times
                        # from 0 to 2
                        splb = True
                    self.reslist[line + 1 + num] = '{:<5s} {:<3} {:>9.6f}   {:>9.6f}   {:>9.6f}   {:>8.4f}     0.04\n' \
                        .format(*resline)
        else:
            for line in id_lines:
                # replacing "REM PART 1 !QNGUYQ2" with "PART 1"
                # restraints should never be placed in this reslist[line]:
                self.reslist[line] = ' '.join(self.reslist[line].split()[1:3]) + '\n'
        # remove PART -1 around AFIX 130:
        for line in neg_part_id_lines:
            self.reslist[line] = ''
        # set refinement cycles back to 8
        shx = ShelxlRefine(self.reslist, self.basefilename, self.fa, self.options)
        shx.set_refinement_cycles('8')
        self.rl.write_resfile(self.reslist, '.res')
        return fatoms

    def cf9(self):
        """
        create disorderd CF3 group on three positions.
        """
        atomline = self.prepare_cf3()
        target_atom = self.dsrp.target[0]
        print('Generating threefold disordered CF3-Group at {}.'.format(target_atom))
        if self.dsrp.dfix:
            restr = dfixr_cf9
        else:
            restr = sadir_cf9
        # make a copy to find fluorine positions:
        ####################################################
        reslist_copy = self.reslist[:]
        self.make_pivot_isotropic(atomline)
        # Fluorine atoms are generated here:
        fatoms = self.make_afix(afixnum=130, linenumber=atomline, resioff=True)
        self.do_refine_cycle(self.rl, self.reslist)
        # this is the bond around the CF3 group rotates
        Y, Z = self.lf.get_bondvector(target_atom)
        Y = remove_partsymbol(Y)
        Z = remove_partsymbol(Z)
        numberedatoms = ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9']
        resistr = ''
        resi0 = ''
        if self.dsrp.resiflag:
            F1, F2, F3, F4, F5, F6, F7, F8, F9 = numberedatoms
            resiclass = self.resi.get_residue_class
            resinum = self.resi.get_resinumber
            resistr = 'RESI ' + resiclass + ' ' + resinum
            resi0 = 'RESI 0\n'
        else:
            nums = NumberScheme(self.reslist, numberedatoms, self.dsrp)
            F1, F2, F3, F4, F5, F6, F7, F8, F9 = nums.get_fragment_number_scheme()
        start_f_coord = self.lf.get_single_coordinate(fatoms[0] + '^A')
        replacelist = (('Z', Z), ('Y', Y), ('F1', F1), ('F2', F2), ('F3', F3),
                       ('F4', F4), ('F5', F5), ('F6', F6), ('F7', F7),
                       ('F8', F8), ('F9', F9))
        for old, new in replacelist:
            restr = [i.replace(old, new) for i in restr]
        ########################################################
        # add restraints to reslist:
        # all fluorines are set, so we can get back to the original res file:
        self.reslist = reslist_copy
        self.rle._reslist = self.reslist
        # Turned out to be problematic with ShelXle:
        # if self.dsrp.occupancy:
        #    occ = self.dsrp.occupancy
        # else:
        occ = str((self.rle.get_fvar_count() + 1) * 10 + 1 + 20)
        fcount = self.rle.get_fvar_count()
        self.rle.set_free_variables(occ, '0.3')
        atomline = self.fa.get_atom_line_numbers([target_atom])[0]
        self.make_pivot_isotropic(atomline)
        atoms_cf9 = ['PART 1 {1}1',
                     F1 + '   {0}   {4}    11.00000    0.04',
                     F2 + '   {0}   {5}    11.00000    0.04',
                     F3 + '   {0}   {6}    11.00000    0.04',
                     'PART 2 {2}1',
                     F4 + '   {0}   {7}    11.00000    0.04',
                     F5 + '   {0}   {8}    11.00000    0.04',
                     F6 + '   {0}   {9}    11.00000    0.04',
                     'PART 3 {3}1',
                     F7 + '   {0}   {10}    11.00000    0.04',
                     F8 + '   {0}   {11}    11.00000    0.04',
                     F9 + '   {0}   {12}    11.00000    0.04',
                     'PART 0 ',
                     self.endm]
        restr = wrap_headlines(restr, 77)
        restr = ''.join(restr).format(fcount + 1, fcount + 2, fcount + 3)
        # place the restraints:
        self.reslist[atomline] = (self.reslist[atomline] + self.startm + restr)
        at1 = self.lf.get_single_coordinate(Y)
        at2 = self.lf.get_single_coordinate(Z)
        coords = []
        # rotate the fluorine coordinate around at1, at2 to get three triples:
        for delta in [0, 120, 240, 40, 160, 280, 80, 200, 320]:
            coord = self.rotate_atom_around_bond(start_f_coord, at1, at2, delta)
            coords.append('{:>10.6f} {:>10.6f} {:>10.6f}'.format(*coord))
        if self.dsrp.resiflag:
            atoms_cf9 = [resistr] + atoms_cf9 + [resi0]
        # join all together:
        atoms_cf9 = '\n'.join(atoms_cf9).format(self.e2s.elem_2_sfac('F'),
                                                fcount + 1, fcount + 2, fcount + 3, *coords)
        self.reslist[atomline] += atoms_cf9
        # have to do this here, because set_free_variables() works on different reslist:
        # self.reslist[self.rle.find_fvarlines()[0]] = ' \n'.join(fvar)+'\n'
        shx = ShelxlRefine(self.reslist, self.basefilename, self.fa, self.options)
        shx.set_refinement_cycles('8')
        self.rl.write_resfile(self.reslist, '.res')
        return fatoms

    def make_afix(self, afixnum, linenumber, resioff=False):
        """
        create an afix to build a CF3 or CH3 group
        :param resioff: bool
        :param linenumber:  line number in the reslist where fragment is placed
        :param afixnum: afix number
        :type afixnum: int
        """
        resistr = ''
        resi0 = ''
        # Turned out to be problematic with ShelXle:
        # if self.dsrp.occupancy:
        #    occ = self.dsrp.occupancy
        # else:
        numberscheme_130 = []
        numberscheme_120 = []
        occ = str((self.rle.get_fvar_count() + 1) * 10 + 1)
        if int(afixnum) == 120:
            self.dsrp.occupancy = occ
            self.rle.set_free_variables(occ, '0.5')
            num_120 = NumberScheme(self.reslist, ['F1', 'F2', 'F3', 'F4', 'F5', 'F6'], self.dsrp)
            # returns also the atom names if residue is active
            numberscheme_120 = num_120.get_fragment_number_scheme()
        if int(afixnum) == 130:
            num_130 = NumberScheme(self.reslist, ['F1', 'F2', 'F3'], self.dsrp)
            # returns also the atom names if residue is active
            numberscheme_130 = num_130.get_fragment_number_scheme()
        if self.dsrp.resiflag and not resioff:
            resiclass = self.resi.get_residue_class
            resinum = self.resi.get_resinumber
            resistr = '\nRESI ' + resiclass + ' ' + resinum
            resi0 = 'RESI 0\n'
            numberscheme_130 = ['F1', 'F2', 'F3']
            numberscheme_120 = ['F1', 'F2', 'F3', 'F4', 'F5', 'F6']
        sfac = self.e2s.elem_2_sfac('F')
        if int(afixnum) == 130:
            # CF3:
            afix_130 = [resistr,
                        'AFIX {0}\nPART -1 !{4}',  # AFIX 120 or 130
                        # 'REM AFIX made by DSR: {3}',
                        numberscheme_130[0] + ' {1}   0 0 0   11.0  0.04',
                        numberscheme_130[1] + ' {1}   0 0 0   11.0  0.04',
                        numberscheme_130[2] + ' {1}   0 0 0   11.0  0.04',
                        # 'REM end of AFIX by DSR {3}', # insert ID and later change it to the PART usw.
                        'PART 0 !{4}\nAFIX 0',
                        resi0]
            afix_130 = '\n'.join(afix_130)
            afix = afix_130
        elif int(afixnum) == 120:
            # CF6:
            # {4} is a second id behind the PART -1 to get rid of it later
            # We need a negative part here to prevent accidental bonding to the pivot atom
            afix_120 = [resistr,
                        'AFIX {0}\nPART -1 !{4}',
                        '\nREM PART 1 !{3}',
                        numberscheme_120[0] + ' {1}  0 0 0   {2}  0.04',
                        numberscheme_120[1] + ' {1}  0 0 0   {2}  0.04',
                        numberscheme_120[2] + ' {1}  0 0 0   {2}  0.04',
                        'REM PART 2 !{3}',
                        numberscheme_120[3] + ' {1}  0 0 0  -{2}  0.04',
                        numberscheme_120[4] + ' {1}  0 0 0  -{2}  0.04',
                        numberscheme_120[5] + ' {1}  0 0 0  -{2}  0.04',
                        'REM PART 0 !{3}',
                        'PART 0 !{4}\nAFIX 0',
                        resi0,
                        self.endm]
            afix_120 = '\n'.join(afix_120)
            afix = afix_120
        else:
            print('Only AFIX 130 and 120 are implemented yet.')
            return False
        # insert the afix (the rand_id is to recognize this line later):
        self.reslist[linenumber] += afix.format(afixnum, sfac, occ, self.rand_id, self.rand_id_part_neg)
        if str(afixnum) == '120':
            return numberscheme_120
        if str(afixnum) == '130':
            return numberscheme_130

    def do_refine_cycle(self, rl, reslist):
        """
        runs a shelxl cycle with L.S. 0
        :param rl: reslist object
        :param reslist: res file as list
        """
        shx = ShelxlRefine(reslist, self.basefilename, self.fa, self.options)
        acta_lines = shx.remove_acta_card()
        shx.set_refinement_cycles('0')
        rl.write_resfile(reslist, '.ins')
        shx.run_shelxl()
        self.lf = ListFile(self.basefilename)
        lst_file = self.lf.read_lst_file()
        shx.check_refinement_results(lst_file)
        rl = ResList(self.res_file)
        reslist = rl.get_res_list()
        shx.restore_acta_card(acta_lines)
        shx.set_refinement_cycles('8')
        rl.write_resfile(reslist, '.res')

    def rotate_atom_around_bond(self, ratom, at1, at2, delta=10):
        """
        R = T**-1*Rx**-1*Ry**-1*Rz*Ry*Rx*T

        1 translate object to origin
        2 rotate around x in xz plane
        3 rotate around y in z axis
        4 rotate aound z with delta angle
        5 transform back with inverse of 1-3

        v = P2 - P1
        |v| = sqrt(vx**2+vy**2+vz**2)
        u = v/|v| = (a, b, c), |u| = 1
        a = (x2-x1)/|v| b = (x2-y1)/|v| c = ...
        """
        ratom = frac_to_cart(ratom, self.cell)
        at1 = frac_to_cart(at1, self.cell)
        at2 = frac_to_cart(at2, self.cell)

        ratom = mpm.matrix(list(ratom) + [1.0])
        delta = radians(delta)
        x0, y0, z0 = at1
        T = mpm.matrix(((1, 0, 0, -x0),
                        (0, 1, 0, -y0),
                        (0, 0, 1, -z0),
                        (0, 0, 0, 1)))
        T1 = mpm.inverse(T)

        vx = at2[0] - at1[0]
        vy = at2[1] - at1[1]  # P2 - P1
        vz = at2[2] - at1[2]
        vnorm = sqrt(vx ** 2 + vy ** 2 + vz ** 2)
        a, b, c = vx / vnorm, vy / vnorm, vz / vnorm
        d = mpm.sqrt(b ** 2 + c ** 2)

        sina = b / d  # For rotation around alpha
        cosa = c / d  #

        Rxa = mpm.matrix(((1, 0, 0, 0),
                          (0, cosa, sina, 0),
                          (0, -sina, cosa, 0),
                          (0, 0, 0, 1)))
        '''
        Rya = mpm.matrix(((      cosa,          0,      -sina,  0),
                         (         0,          1,          0,  0),
                         (      sina,          0,       cosa,  0),
                         (         0,          0,          0,  1)))
        ''' '''
        Rza = mpm.matrix(((      cosa,       sina,          0,  0),
                         (     -sina,       cosa,          0,  0),
                         (         0,          0,          1,  0),
                         (         0,          0,          0,  1)))
        '''
        Rxa1 = mpm.inverse((Rxa))

        cosb = d  #
        sinb = -a  # rotation around beta

        Ryb = mpm.matrix(((cosb, 0, -sinb, 0),
                          (0, 1, 0, 0),
                          (sinb, 0, cosb, 0),
                          (0, 0, 0, 1)))

        Ryb1 = mpm.inverse((Ryb))

        sind = mpm.sin(delta)
        cosd = mpm.cos(delta)

        Rzd = mpm.matrix(((cosd, sind, 0, 0),
                          (-sind, cosd, 0, 0),
                          (0, 0, 1, 0),
                          (0, 0, 0, 1)))

        R = T1 * Rxa * Ryb * Rzd * Ryb1 * Rxa1 * T
        v = R * ratom
        v = cart_to_frac(v[:3], self.cell)
        v = [round(i, 6) for i in v]
        return v

    def make_cf3_thorus(self, atom=None):
        """
        This method was just an experiment. The CF9 method is sufficient.

        Creates a thorus of isotropic fluorine atoms around the central
        atom of a cf3 group. The occupancy is estimated from the
        residual density values in the lst file.

        :param atom: central atom of the cf3 group
        :type atom: string
        """
        coords = []
        if not atom:
            atom = self.dsrp.target[0]
        # returns the atom names of the fluorine atoms:
        fluorine_names = self.cf3()
        # self.do_refine_cycle(self.rl, self.reslist)
        print(fluorine_names[0], '##################')
        ratom = self.lf.get_single_coordinate(fluorine_names[0])
        print('ratom:', ratom)
        # ratom = self.get_coordinates_of_first_atom(fluorine_names)
        self.lf.read_lst_file()
        bondvec = self.lf.get_bondvector()
        print('bondvec:', bondvec)
        at1 = self.lf.get_single_coordinate(bondvec[0].upper())
        at2 = self.lf.get_single_coordinate(bondvec[1].upper())
        names = ['F{}'.format(i) for i in range(1, 25)]
        num_thor = NumberScheme(self.reslist, names, self.dsrp)
        names = num_thor.get_fragment_number_scheme()
        for delta in range(0, 360, 15):
            coord = self.rotate_atom_around_bond(ratom, at1, at2, delta)
            coords.append(coord)
        # print('PART 0')
        diffden = self.lf.get_difference_density(averaged=False)
        print(diffden)
        rotation = self.lf.get_degree_of_highest_peak()
        print(rotation)
        n = (rotation / 15) + 4  # +1
        diffden = shift(diffden, n)
        print(len(diffden))
        print('\n')
        print('')
        print('AFIX 6')
        # make restraints:
        sad12 = []
        sad13 = []
        ff = []
        flast = False
        for i, num, co, dif in zip(names, [i for i in range(1, 25)], coords, diffden):
            # print('PART {}'.format(num))
            if not flast:
                flast = 'F24'
            ff.append((flast, i))
            flast = i
            sad12.append(('C22', i))
            sad13.append(('C19', i))
            print('{}  3  {:0<8.6}  {:0<8.6}  {:0<8.6}  {:<8.6}  -1.2\
            '.format(i, co[0], co[1], co[2], 10.0 + (dif / 2760.0)))
            # with free variables:
            # print('{}  3  {:0<8.6}  {:0<8.6}  {:0<8.6}  {:<8.6}  261\
            # '.format(i, co[0], co[1], co[2], 10.0*num+1+10))
        from misc import flatten
        print('SADI ' + ' '.join(flatten(sad12)))
        print('SADI ' + ' '.join(flatten(sad13)))
        print('SADI ' + ' '.join(flatten(ff)))
        d = 0
        for dif in diffden:
            d += dif / 2760.0
        # print('PART 0')
        print('AFIX 0')
        print(d)


if __name__ == '__main__':
    pass
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


def calc_ellipsoid_axes(coords, uvals, cell, probability=0.5, longest=True):
    """
    This method calculates the principal axes of an ellipsoid as list of two
    fractional coordinate triples.
    Many thanks to R. W. Grosse-Kunstleve and P. D. Adams
    for their great publication on the handling of atomic anisotropic displacement
    parameters:
    R. W. Grosse-Kunstleve, P. D. Adams, J Appl Crystallogr 2002, 35, 477–480.

    F = ... * exp ( -2π²[ h²(a*)²U11 + k²(b*)²U22 + ... + 2hka*b*U12 ] )

    SHELXL atom:
    Name type  x      y      z    occ     U11 U22 U33 U23 U13 U12
    F3    4    0.210835   0.104067   0.437922  21.00000   0.07243   0.03058 =
       0.03216  -0.01057  -0.01708   0.03014
    >>> import src.mpmath as mpm
    >>> cell = [10.5086, 20.9035, 20.5072, 90, 94.13, 90]
    >>> coords = [0.210835,   0.104067,   0.437922]
    >>> uvals = [0.07243, 0.03058, 0.03216, -0.01057, -0.01708, 0.03014]
    >>> l = calc_ellipsoid_axes(coords, uvals, cell, longest=True)
    >>> print(mpm.nstr(l))
    [[0.24765096, 0.11383281, 0.43064756], [0.17401904, 0.09430119, 0.44519644]]
    >>> calc_ellipsoid_axes(coords, uvals, cell, longest=False)
    [[[0.24765096, 0.11383281, 0.43064756], [0.218406, 0.09626142, 0.43746127], [0.21924358, 0.10514684, 0.44886868]], [[0.17401904, 0.09430119, 0.44519644], [0.203264, 0.11187258, 0.43838273], [0.20242642, 0.10298716, 0.42697532]]]
    >>> cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
    >>> coords = [0.210835,   0.104067,   0.437922]
    >>> uvals = [0.07243, -0.03058, 0.03216, -0.01057, -0.01708, 0.03014]
    >>> calc_ellipsoid_axes(coords, uvals, cell, longest=True)
    <BLANKLINE>
    Ellipsoid is non positive definite!
    <BLANKLINE>
    False

    >>> uvals = [0.07243, 0.03058, 0.03216, -0.01057, -0.01708]
    >>> calc_ellipsoid_axes(coords, uvals, cell, longest=False)
    Traceback (most recent call last):
    ...
    Exception: 6 Uij values have to be supplied!

    >>> cell = (10.5086, 20.9035, 90, 94.13, 90)
    >>> coords = [0.210835,   0.104067,   0.437922]
    >>> uvals = [0.07243, 0.03058, 0.03216, -0.01057, -0.01708, 0.03014]
    >>> calc_ellipsoid_axes(coords, uvals, cell, longest=True)
    Traceback (most recent call last):
    ...
    Exception: cell needs six parameters!

    :param coords: coordinates of the respective atom in fractional coordinates
    :type coords: list
    :param uvals: Uij valiues of the respective ellipsoid on fractional
                  basis like in cif and SHELXL format
    :type uvals: list
    :param cell: unit cell of the structure: a, b, c, alpha, beta, gamma
    :type cell:  list
    :param probability: thermal probability of the ellipsoid
    :type probability: float or int
    :param longest: not always the length is important. make to False to
                    get all three coordiantes of the ellipsoid axes.
    :type longest: boolean

    """
    probability += 1
    # Uij is symmetric:
    if len(uvals) != 6:
        raise ValueError('6 Uij values have to be supplied!')
    if len(cell) != 6:
        raise ValueError('cell needs six parameters!')
    # orthogonalization matrix that transforms the fractional coordinates
    # with respect to a crystallographic basis system to coordinates
    # with respect to a Cartesian basis:
    _A = A(cell).orthogonal_matrix
    Ucart = ufrac_to_ucart(_A, cell, uvals)
    # print(Ucart)
    # E => eigenvalues, Q => eigenvectors:
    E, Q = mpm.eig(mpm.matrix(Ucart))
    # calculate vectors of ellipsoid axes
    try:
        sqrt(E[0])
        sqrt(E[1])
        sqrt(E[2])
    except ValueError:
        print('\nEllipsoid is non positive definite!\n')
        return False
    v1 = mpm.matrix([Q[0, 0], Q[1, 0], Q[2, 0]])
    v2 = mpm.matrix([Q[0, 1], Q[1, 1], Q[2, 1]])
    v3 = mpm.matrix([Q[0, 2], Q[1, 2], Q[2, 2]])
    v1i = v1 * (-1)
    v2i = v2 * (-1)
    v3i = v3 * (-1)
    # multiply probability (usually 50%)
    e1 = sqrt(E[0]) * probability
    e2 = sqrt(E[1]) * probability
    e3 = sqrt(E[2]) * probability
    # scale axis vectors to eigenvalues
    v1, v2, v3, v1i, v2i, v3i = v1 * e1, v2 * e2, v3 * e3, v1i * e1, v2i * e2, v3i * e3
    # find out which vector is the longest:
    length = mpm.norm(v1)
    v = 0
    if mpm.norm(v2) > length:
        length = mpm.norm(v2)
        v = 1
    elif mpm.norm(v3) > length:
        length = mpm.norm(v3)
        v = 2
    # move vectors back to atomic position
    atom = mpm.matrix(_A) * mpm.matrix(coords)
    v1, v1i = v1 + atom, v1i + atom
    v2, v2i = v2 + atom, v2i + atom
    v3, v3i = v3 + atom, v3i + atom
    # go back into fractional coordinates:
    a1 = cart_to_frac(v1, cell)
    a2 = cart_to_frac(v2, cell)
    a3 = cart_to_frac(v3, cell)
    a1i = cart_to_frac(v1i, cell)
    a2i = cart_to_frac(v2i, cell)
    a3i = cart_to_frac(v3i, cell)
    allvec = [[a1, a2, a3], [a1i, a2i, a3i]]
    if longest:
        # only the longest vector
        return [allvec[0][v], allvec[1][v]]
    else:
        # all vectors:
        return allvec
