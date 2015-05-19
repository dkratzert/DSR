#-*- encoding: utf-8 -*-
#möp
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <daniel.kratzert@ac.uni-freiburg.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
# each residue number must always have the same class!
# e.g 1 can be only TOL, not also NET4 but Tol can be 1 or 2.
#
# Wenn das dict leer ist werden überhaupt keine RESI betreffenden dinge in das insfile übernommen

from __future__ import print_function
import sys
from resfile import ResList
from constants import RESTRAINT_CARDS


class Resi(object):
    def __init__(self, reslist, dsr_line_dict, dbhead, db_residue_string, find_atoms):
        '''
        Handles the RESI instructions and restraints
        self._dsr_command_resi_list : list of RESI commands from command line.
        self._resi_dict_dsr_command : dictionary of commands after get_resi_syntax()

        self._dsr_dict['resi'] = False : no residue given
        self._dsr_dict['resi'] = '' : only resi command given
        self._dsr_dict['resi'] = ['number, e.g. 3'] : only resi number given
        self._dsr_dict['resi'] = ['class'] : only resi class given
        self._dsr_dict['resi'] = ['class', 'number'] : resi class and number given
        self._dsr_dict['resi'] = ['class', 'number', 'alias'] : resi class, number and alias given
        :param reslist: res file as list
        :type reslist: list
        :param dsr_line_dict: dsr command line as dictionary
        :type dsr_line_dict: dictionary
        :param dbhead: database header with restraints and residue
        :type dbhead: list
        :param db_residue_string: db residue string
        :type db_residue_string: string
        :param find_atoms: find_atoms object
        :type find_atoms: object
        '''
        self._reslist = reslist
        self._find_atoms = find_atoms
        self._dbhead = dbhead
        self._atoms_in_reslist = self._find_atoms.collect_residues()
        self._residues_in_res = sorted(self._atoms_in_reslist.keys())
        self._dsr_dict = dsr_line_dict.copy()
        self._dsr_command_resi_list = self._dsr_dict['resi']
        if self._dsr_command_resi_list:
            # at least one residue parameter defined in .res file
            self._resi_dict_dsr_command = self.get_resi_syntax(self._dsr_command_resi_list)
        if self._dsr_command_resi_list == '':
            # use residue from database
            self._resi_dict_dsr_command = {'class': None, 'number': None, 'alias': None}
        if self._dsr_command_resi_list == False:
            # residues tuned off
            self._resi_dict_dsr_command = False
        try:
            # use residue class from db if not given in dsr command:
            self._db_resi_list = db_residue_string.split()
        except(AttributeError):
            print('No valid residue "RESI classname" in the database entry '\
                    'of {} found.'.format(self._dsr_dict['fragment']))
            sys.exit()
        self._resi_dict_db = self.get_resi_syntax(self._db_resi_list)
        self._combined_resi = self.build_up_residue()

    @property
    def get_resinumber(self):
        '''
        Returns the residue number of the currently fitted fragment
        :type self._combined_resi['number']: string
        '''
        return self._combined_resi['number']

    @property
    def get_residue_class(self):
        '''
        Returns the residue class of the currently fitted fragment
        :type self._combined_resi['class']: string
        '''
        return self._combined_resi['class']

    def remove_resi(self, head):
        '''
        removes all resi commands and classes from head
        :param head: database header of the current fragment
        :type head: list
        '''
        rhead = [] #head without resi
        delhead = []
        for dummy, line in enumerate(head):
            line = line.split()
            try:
                if line[0][:4] in RESTRAINT_CARDS:
                    line[0] = line[0].split('_')[0]
            except(IndexError):
                continue
            line = ' '.join(line)
            delhead.append(line)
        for line in delhead:
            line = line.strip(' \n\r').upper()
            if line.startswith('RESI'):
                    continue
            rhead.append(line)
        return rhead

    def format_restraints(self, head):
        '''
        in case of RESI, format the restraints like "SAME_class"
        '''
        newhead = []
        for line in head:
            line = line.upper()
            line = line.split()
            try:
                line[0]
            except:
                continue
            if line[0] in RESTRAINT_CARDS:
                line[0] = line[0]+'_'+self._combined_resi['class']
                line = ' '.join(line)
            else:
                line = ' '.join(line)
            newhead.append(line)
        return newhead


    def get_unique_resinumber(self, resinum):
        '''
        Finds a unique resi number. If the number is already unique
        the given is used.
        :param resinum: residue number of the fragment
        :type resinum: string
        :return resinum: unique residue number
        '''
        new_num = '1'
        if not resinum:
            while new_num in self._residues_in_res:
                new_num = str(int(new_num)+1)
            return new_num
        if resinum in self._residues_in_res:
            print('Warning: The residue number "{}" you have chosen is already '\
                    'in use!'.format(resinum))
            while new_num in self._residues_in_res:
                new_num = str(int(new_num)+1)
            return new_num
        else:
            return resinum


    def build_up_residue(self):
        '''
        Decides which class and residue number should be used for the fragment.
        Returns a final dict with the residue settings.
        self._resi_dict_dsr_command is False if resi is enabled but no values given.
        '''
        final_residue = {'class': None, 'number': None, 'alias': None}
        resiclass = None
        resinum = None
        resialias = None

        #### for the db entry
        if self._resi_dict_db['alias']:
            resialias = self._resi_dict_db['alias']
        if self._resi_dict_db['class']:
            resiclass = self._resi_dict_db['class'].upper()
        if self._resi_dict_db['number']:
            resinum = self._resi_dict_db['number']
        #### for the comlist entry
        if self._resi_dict_dsr_command == False:
            return final_residue
        else:
            if self._resi_dict_dsr_command['alias']:
                resialias = self._resi_dict_dsr_command['alias']
            if self._resi_dict_dsr_command['class']:
                resiclass = self._resi_dict_dsr_command['class'].upper()
            if self._resi_dict_dsr_command['number']:
                resinum = self._resi_dict_dsr_command['number']
        final_residue = {'class': resiclass, 'number': resinum, 'alias': resialias}
        if not resinum:
            final_residue['number'] = self.get_unique_resinumber(resinum)
            print('No residue number was given. Using residue number {}.'.format(final_residue['number']))
        else:
            final_residue['number'] = self.get_unique_resinumber(resinum)
            print('Using residue number {}.'.format(final_residue['number']))
        if final_residue['class']:
            return final_residue


    def _wrong_syntax(self):
        print('This is not a valid RESIdue syntax!')


    def get_resi_syntax(self, resi):
        '''
        Checks if resi class, number and/or alias are present and valid.
        start with digit-> rest auch digit-> resinumber or alias
        start with letter-> rest letter or digit -> residue class

        The return value of this method is {'class', 'number', 'alias'}
        as a dictionary. Alias is empty if not given.

        The return value of just "RESI" in the command line is an empty dict
        :param resi: residue definition like ['3', 'CF3']
        :type resi: list
        '''
        resi_dict = {
            'class' : None,
            'number': None,
            'alias' : None}

        if not resi:
            print('No valid RESI instruction found in the database entry!')
            sys.exit()

        try:
            resi.sort()
        except(AttributeError):
            return resi_dict

        if str.isalpha(resi[-1][0]): # first character of class must be a letter
            resi_dict['class'] = resi.pop()
        if len(resi) > 0:
            if str.isdigit(resi[0][0]): # first character of number must be an digit
                resi_dict['number'] = resi[0]
                del resi[0]
            else:
                self._wrong_syntax()
                # not totally correct, but less difficult:
                print('Only digits are allowed for the residue number.')
                sys.exit()
        if len(resi) > 0:
            if str.isdigit(resi[0]):
                resi_dict['alias'] = resi[0]
                del resi[0]
            else:
                self._wrong_syntax()
                print('Only digits are allowed in the residue alias.')
                sys.exit()
        if resi_dict['number']:
            if len(resi_dict['number']) > 4:
                self._wrong_syntax()
                print('Only four digits allowed in residue number!')
                sys.exit()
        return resi_dict


    def resi_class_atoms_consistent(self):
        '''
        currently not used

        find out if the atom names of the current residue class fit
        to the atom names in already existing classes.
        
        go through residues:
            go through atoms:
                if atom == current residue class:
                    add atom to at_list
            compare for residue if fragment atom list fits to at_list
        '''
        for num in list(self._atoms_in_reslist.keys()):
            print(num, len(self._atoms_in_reslist[num]), self._atoms_in_reslist[num][:][0][3], \
                    [i[0] for i in self._atoms_in_reslist[num][:]])





if __name__ == '__main__':
    from dsrparse import DSR_Parser
    from dbfile import global_DB
    from atomhandling import FindAtoms
    from resfile import ResListEdit
    res_file = 'p21c.res'
    #dbhead = ['REM test\n', 'RESI 1 TOL\n', 'SADI_TOL C1 C2\n']
    rl = ResList(res_file)
    res_list = rl.get_res_list()
    find_atoms = FindAtoms(res_list)
    rle = ResListEdit(res_list, find_atoms)
    dsrp = DSR_Parser(res_list, rle)
    dsr_dict = dsrp.get_dsr_dict
    #fragment = dsr_dict['fragment']
    fragment = 'toluene'
    invert = True
    gdb = global_DB(invert)
    db = gdb.build_db_dict()
    fragline = gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
    dbatoms = gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
    dbhead = gdb.get_head_from_fragment(fragment)        # this is only executed once

    for i in dbhead:
        print(i)
    residue = '5 CCF3'
    resi = Resi(res_list, dsr_dict, dbhead, residue, find_atoms)
    # db_resi = resi.get_resi_from_db()
    # print 'the new resinumber:', resi.get_unique_resinumber()
    print()
    print()
    head = resi.make_resihead()
    #for i in head:
    #    print(i)

    resiatoms = find_atoms.collect_residues()
    print(resiatoms['1'])
    #if resi:
    #    print 'ja, resi aktivieren! hallo'
    #else:
    #    print 'nein, resi nicht aktivieren!'
