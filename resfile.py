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
from __future__ import print_function

import textwrap

import misc
import sys
import os


__metaclass__ = type  # use new-style classes


def filename_wo_ending(resfilename):
    '''
    returns the input file name without ending
    :param resfilename: string like 'p21c.res'
    '''
    if not resfilename:
        return False
    file_ext = os.path.splitext(resfilename)
    basefile = file_ext[0]
    if not file_ext[1] in ['.res', '']:
        print("Please give a .res file as argument.")
        sys.exit(0)
    return str(basefile)


class ResList():
    '''
    Reads and writes the res-file as list
    '''
    def __init__(self, res_file_name):
        '''
        :param res_file_name: string like 'p21c.res'
        '''
        self.__res_file_name = res_file_name
        self.__basefile = filename_wo_ending(self.__res_file_name)


    def get_res_list(self):
        '''
        read the res-file and return it as list
        '''
        filename = self.__basefile+'.res'
        reslist = []
        try:
            with open(filename, 'r') as f:
                for line in f:
                    reslist.append(line)
        except(IOError):
            print('Unable to read {} file.'.format(filename))
            sys.exit()
        return reslist


    def write_resfile(self, reslist, ending):
        '''
        Writes the content of reslist to basefile+ending.
        Ending should be '.ins' or '.res'
        :param reslist: SHELXL .res file as list
        :param ending:  string, file ending like '.res'
        '''
        try:
            nfile = open(self.__basefile+ending, 'w')  # open the ins file
        except(IOError):
            print('Unable to write .res file! Check directory write permissions.')
            sys.exit(-1)
        for i in reslist:            #modified reslist
            if i == ' ':
                i = i+' \n'  # prevent space character coming in front of an atom accidentally
            nfile.write("%s" %i)    #write the new file
        nfile.close()


class ResListEdit():
    '''
    Class to add, remove and list lines of the resfile
    '''
    def __init__(self, reslist, find_atoms):
        '''

        :param reslist: SHELXL .res file as list
        :param find_atoms: FindAtoms() object
        '''
        self._reslist = reslist
        self._find_atoms = find_atoms
        self.fvarlines = self.find_fvarlines()
    
    def get_cell(self):
        '''
        returns cell parameters of the res file
        '''
        cell = []
        for line in self._reslist:
            if line.startswith('CELL'):
                cell = line.split()[2:8]
                cell = [float(i) for i in cell]
        return cell

    def add_line(self, linenum, new_line):
        '''
        adds a single line of string to the res file
        checks are missing which make sure the data structure of new_line is correct!
        :param linenum:  integer, line number to add string
        :param new_line: new string to insert like 'foo bar'
        '''
        linenum = linenum-1
        self._reslist.insert(linenum, new_line+'\n')
        return self._reslist


    def remove_line(self, linenum, rem=False, remove=False, frontspace=False):
        '''
        removes a single line from the res file with tree different methods.
        The default is a space character in front of the line (frontspace).
        This removes the line in the next refinement cycle. "rem" writes rem
        in front of the line and "remove" clears the line.
        :param linenum: integer, line number
        :param rem:     True/False, activate comment with 'REM' in front
        :param remove:  True/False, remove the line
        :param frontspace: True/False, activate removing with a front space
        '''
        line = self._reslist[linenum]
        if rem:   # comment out with 'rem ' in front
            self._reslist[linenum] = 'rem '+line
            if misc.multiline_test(line):
                self._reslist[linenum+1] = 'rem '+self._reslist[linenum+1]
        elif remove:  # really delete the line "linenum"
            if misc.multiline_test(line):
                self._reslist[linenum] = ''
                self._reslist[linenum+1] = ''
            else:
                self._reslist[linenum] = ''
        if frontspace:  # only put a space in front
            self._reslist[linenum] = ' '+line
            if misc.multiline_test(line):
                self._reslist[linenum+1] = ' '+self._reslist[linenum+1]
        return self._reslist


    def list_lines(self, startline, endline):
        '''
        returns the lines between startline and endline as list
        :param startline: integer
        :param endline:   integer
        '''
        lines = []
        try:
            start = int(startline)
            end = int(endline)
        except(NameError, SyntaxError):
            return False
        for line, i in enumerate(self._reslist):
            if line >= start:
                lines.append(i)
            if line == end:
                break
        return lines


    def getAll(self):
        '''returns the whole resfile as list'''
        return self._reslist


    def find_fvarlines(self):
        '''
        Finds the FVAR line or the first line with an atom.
        returns a list of FVAR occurences
        sys.exit() if no atom or FVAR found
        '''
        fvarlines = misc.find_multi_lines(self._reslist, r'^FVAR')
        if not fvarlines:   # There is no FVAR in the res file after SHELXS!
            for num, i in enumerate(self._reslist):
                if self._find_atoms.is_atom(i):
                    first_atom = num
                    break
            if not first_atom:
                print('\nNo atom or Q-peak found! Can not proceed...\n')
                sys.exit()
            fvarlines = []
            fvarlines.append(first_atom)
            self._reslist.insert(first_atom, ' \n')
            return fvarlines
        else:
            return fvarlines

    def insert_frag_fend_entry(self, dbatoms, fragline, fvarlines):
        '''
        Inserts the FRAG ... FEND entry in the res file.
        :param dbatoms:   list of atoms in the database entry
        :param fragline:  string with "FRAG 17 cell" from the database entry
        :param fvarlines: line where FVAR or the first atom is located
        '''
        dblines = []
        db = [list(map(str, i)) for i in dbatoms]
        for i in db:
            i[2] = '{:>10.6f}'.format(float(i[2]))
            i[3] = '{:>10.6f}'.format(float(i[3]))
            i[4] = '{:>10.6f}'.format(float(i[4]))
            dblines.append('    '.join(i).rstrip())
        dblines = '\n'.join(dblines)
        dblines = '  '.join(fragline)+'\n'+dblines
        dblines = '\n The following is from DSR:\n'+dblines
        dblines = dblines+'\nFEND\n\n'
        #dblines = misc.ll_to_string(db)+'\n'
        self._reslist.insert(fvarlines[-1]+1, dblines)   # insert the db entry right after FVAR

    def get_fvarlist(self):
        fvar_list = []
        for line in self.fvarlines:
            fvar = self._reslist[line].split()
            if fvar:
                del fvar[0]
            fvar_list.extend(fvar)
        return fvar_list

    def set_free_variables(self, occupancynumber, fvalue='0.5'):
        '''
        Inserts additional free variables according to the occ parameter
        This function starts at the end of parse_dsr_line() so we don't have
        to care about it anywhere else.
        :param occupancynumber: string, like '21.0'
        :type occupancynumber: string
        :param fvarlines:       list, list of line numbers where FVAR is located
                                      in the res file
        '''
        fvar_list = self.get_fvarlist()  # free variables in the res file
        varlen = self.get_fvar_count()
        occupancynumber = occupancynumber.strip('-')
        # how many numbers do we have?:
        # the occupancynumber is split in the fvar part and the occupancy part:
        num = occupancynumber.split('.')   
        fvar = int(num[0])//10       # e.g. 20.5 is fvar 2 and occupancy 0.5
        if fvar == 0:
            fvar = 1
        difference = (fvar - varlen)
        if fvar > 1:
            if difference > 0:
                for i in range(difference):
                    fvar_list.append(fvalue)   # if an fvar is missing, add a new one
            else: # make sure the occupancy of the disoerder parts get not < 0:
                if len(fvar_list)-(fvar-1) >= 0: # make sure
                    fvar_value = fvar_list[fvar-1]
                    if (float(occupancynumber)-(10*int(fvar)+float(fvar_value))) < 0:
                        fvar_list[fvar - 1] = 0.5
        lines = misc.chunks(fvar_list, 7)
        if len(fvar_list) != 0:
            for line in self.fvarlines:
                self.remove_line(line, remove=True) # removes the old FVAR
        fvars = [' '.join(i) for i in lines]
        fvars = ['FVAR '+i for i in fvars]
        self._reslist[self.fvarlines[0]] = ' \n'.join(fvars)
        return fvars

    def get_fvar_count(self):
        '''
        returns the last used free variable defined with FVAR
        '''
        fvars = len(self.get_fvarlist())
        if not fvars:
            fvars = 0
        return fvars

# for testing
if __name__ == '__main__':
    from dbfile import global_DB
    from atomhandling import FindAtoms
    invert = True
    res_file = 'p21c.res'
    res_list = ResList(res_file)
    reslist =  res_list.get_res_list()
    find_atoms = FindAtoms(reslist)
    rle = ResListEdit(reslist, find_atoms)
    fvarlines = rle.find_fvarlines()
    gdb = global_DB(invert)
    db = gdb.build_db_dict()
    fragment = 'toluene'
    dbatoms = gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
    dbhead = gdb.get_head_from_fragment(fragment)        # this is only executed once

    fragline = gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
    rle.insert_frag_fend_entry(dbatoms, fragline, fvarlines)
    rle.set_free_variables('91.1234')

    for num, i in  enumerate(reslist):
        print(num+1, ''.join(i.strip('\n\r')))
        if num > 15:
            break
    print()

    reslist = rle.add_line(3, 'Hallo Welt!!!!!!!!!!!!!')

    print()
    for num, i in  enumerate(reslist):
        print(num+1, ''.join(i.strip('\n')))
        if num > 15:
            break


    reslist = rle.add_line(2, 'xghrzdg###################fdghgfdgh')


    print()
    for n, i in enumerate(reslist):
        print(n+1, ''.join(i.strip('\n')))
        if n > 5:
            break

    res_list.write_resfile(reslist, '.tst')

    print('\n########################################\n')

    for i in rle.list_lines(15, 33):
        print(i.strip('\n'))

    print(misc.ll_to_string(dbatoms))


