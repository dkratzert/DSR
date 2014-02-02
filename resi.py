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
    '''Handles the RESI instructions and restraints
    reslist: the resfile as list
    dsr_line: the resi command from the dsr_line. e.g. RESI TOL or RESI 2 TOL
    dbhead: db header from dbfile module
    
    self.__com_resi_list : list of RESI commands from command line. 
                          self.__com_resi_list is == 'DB' if only RESI is given in res file.
    self.__resi_dict_com : dictionary of commands after get_resi_syntax()
                If self.__com_resi_list is == 'DB' then self.__resi_dict_com == False
    '''

    def __init__(self, reslist, dsr_line, dbhead, residue, find_atoms):
        self.__reslist = reslist
        self._find_atoms = find_atoms
        self._atoms_in_reslist = self._find_atoms.collect_residues()
        self._residues_in_res = sorted(self._atoms_in_reslist.keys())
        self.__dsr_dict = dsr_line.copy()
        self.__command = self.__dsr_dict['command']
        self.__com_resi_list = self.__dsr_dict['resi'][:] # makes a copy because we need it also later
        if self.__com_resi_list:
            if self.__com_resi_list != 'DB':
                ########################################
                self.__resi_dict_com = self.get_resi_syntax(self.__com_resi_list)
                #########################################
        else: 
            self.__resi_dict_com = False
        if self.__com_resi_list == 'DB':
            self.__resi_dict_com = 'Empty'
        self.__dbhead = dbhead
        #self.__db_resi_list = self.get_resi_from_db() #new header has nor residue, but db_dict has resi key
        self.__db_resi_list = residue.split()
        ##############################
        self.__resi_dict_db = self.get_resi_syntax(self.__db_resi_list)
        ##############################
        self.__combined_resi = self.build_up_residue()
        #print '\nFinales dict fuer resi:', self.__combined_resi


    def remove_resi(self, head):
        '''removes all resi commands and classes from head'''
        rhead = [] #head without resi
        delhead = []
        for num, line in enumerate(head):
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



    def make_resihead(self):
        '''
        Finally builds up the afix header with residue informations.
        
        '''
        resi_to_insert = []
        resi = False
        #we need at least a class
        if self.__dsr_dict['resi'] or self.__combined_resi['class']:
            resi = True
        if resi: 
            head = self.remove_resi(self.__dbhead)
            residef = (['RESI', self.get_unique_resinumber(), 
                        self.__combined_resi['class'], 
                        self.__combined_resi['alias']])
            for i in residef:
                if i:
                    resi_to_insert.append(i)
            residef = ' '.join(resi_to_insert)
            head.append(residef)
            head = self.format_restraints(head)
            return head
        else:
            # residues deactivated, removing resi from head
            head = self.remove_resi(self.__dbhead)
            return head
    
    
    
    def format_restraints(self, head):
        '''
        in case of RESI, format the restraints like "SAME_class"
        '''
        newhead = []
        for line in head:
            line = line.upper()
            line = line.split()
            if line[0]in RESTRAINT_CARDS:
                line[0] = line[0]+'_'+self.__combined_resi['class']
                line = ' '.join(line)
            else:
                line = ' '.join(line)
            newhead.append(line)
        #print newhead
        return newhead
    
            
    def get_unique_resinumber(self):
        '''
        Finds a unique resi number. If the number is already unique
        the given is used.
        '''
        try:
            resinum = self.__combined_resi['number']
        except(AttributeError):
            resinum = False
        new_num = '1'
        if not resinum:
            while new_num in self._residues_in_res:
                new_num = str(int(new_num)+1)
            return new_num
        if resinum in self._residues_in_res:
            print('Warning: The residue number "{}" you have chosen is already in use!'.format(resinum))
            while new_num in self._residues_in_res:
                new_num = str(int(new_num)+1)
            print('         I am using number "{}" instead.\n'.format(new_num))
            return new_num
        else:
            return self.__combined_resi['number']    
    
    
    def build_up_residue(self):
        '''
        Decides which class and resideue number should be used for the fragment. 
        Returns a final dict with the residue settings.
        '''
        final_residue = {'class': None, 'number': None, 'alias': None}
        resiclass = None
        resinum = None
        resialias = None
        
        #### for the db entry
        if self.__resi_dict_db['alias']:
            resialias = self.__resi_dict_db['alias']
        if self.__resi_dict_db['class']:
            resiclass = self.__resi_dict_db['class']
        if self.__resi_dict_db['number']:
            resinum = self.__resi_dict_db['number']
        #### for the comlist entry    
        if not self.__resi_dict_com:
            return final_residue
        if self.__resi_dict_com != 'Empty':
            if self.__resi_dict_com['alias']:
                resialias = self.__resi_dict_com['alias']
            if self.__resi_dict_com['class']:
                resiclass = self.__resi_dict_com['class']
            if self.__resi_dict_com['number']:
                resinum = self.__resi_dict_com['number']
        final_residue = {'class': resiclass, 'number': resinum, 'alias': resialias}
        if final_residue['class']:# and final_residue['number']:
            return final_residue
    
 
    
    def get_resi_from_db(self):
        '''
        gets the RESI num class from the db
        '''
        resi_list = False
        for line in self.__dbhead:
            line = line.upper()
            if line.startswith('RESI'):
                resi_list = line.split()
                del resi_list[0]
                return resi_list
                break
        if not resi_list:
            print('No valid RESI instruction found in the database entry.')
            sys.exit()
        
    
    
    def _wrong_syntax(self):
        print('This is not a valid RESIDUE syntax!')


    def get_resi_syntax(self, resi):
        '''
        Checks if resi class, number and/or alias are present and valid.
        start with digit-> rest auch digit-> resinumber or alias
        start with letter-> rest letter or digit -> residue class

        The return value of this method is {'class', 'number', 'alias'}
        as a dictionary. Alias is empty if not given. 
        
        Is it important which one is the number and which one the alias?
        
        The return value of just "RESI" in the command line is an empty dict
        '''
        resi_dict = {
            'class' : None, 
            'number': None, 
            'alias' : None}
        
        if resi == None:
            print('No valid RESI instruction found in the database entry!')
            sys.exit()
        else:
            resi.sort()
            if str.isalpha(resi[-1][0]):
                resi_dict['class'] = resi.pop()
            if len(resi) > 0:
                if str.isdigit(resi[0]):
                    resi_dict['number'] = resi[0]
                    del resi[0]
                else:
                    self._wrong_syntax()
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
            else: 
                if self.__com_resi_list == 'DB': # in this case no residue number is given at all.
                    number = self.get_unique_resinumber()
                    #number = 3
                    print('No residue number was given. Using residue number {}.'.format(number))
                #self.get_unique_resinumber()
                #resi_dict['number'] = True
            return resi_dict


    def resi_class_atoms_consistent(self):
        '''
        find out if the atom names of the current residue class fit
        to the atom names in already existing classes.
        '''
        pass
        

    def get_unique_residue_name(self):
        '''if not given, finds a unique RESIdue name'''
        pass


    def multi_residues_single_fragment(self):
        '''if exist RESI 1 Fragment:
                add RESI 2 Fragment
        '''
        pass





if __name__ == '__main__':
    from options import OptionsParser
    from dbfile import global_DB
    options = OptionsParser()
    #dbhead = ['REM test\n', 'RESI 1 TOL\n', 'SADI_TOL C1 C2\n']
    rl = ResList(options.res_file)
    res_list = rl.get_res_list()
    dsrp = DSR_Parser(res_list)
    dsr_dict = dsrp.parse_dsr_line()
    fragment = dsr_dict['fragment']
    gdb = global_DB()
    db = gdb.build_db_dict()
    fragline = gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
    dbatoms = gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
    dbhead = gdb.get_head_from_fragment(fragment)        # this is only executed once

    
    resi = Resi(res_list, dsr_dict, dbhead)
    db_resi = resi.get_resi_from_db()
   # print 'the new resinumber:', resi.get_unique_resinumber()
    print()
    print()
    print(resi.make_resihead())
    
    #if resi:
    #    print 'ja, resi aktivieren! hallo'
    #else:
    #    print 'nein, resi nicht aktivieren!'
