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
import misc
import shutil
import re
import sys
import os
import constants


__metaclass__ = type  # use new-style classes

def get_cell(res_list):
    '''
    Returns the unit cell parameters from the res file as list:
    ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
    '''
    cell = False
    for line in res_list:
        if line.startswith('CELL'):
            cell = line.split()[2:]
            break
    if not cell:
        print('Unable to find unit cell parameters in th res file.')
        sys.exit()
    return cell

class ResList():
    '''Reads and writes the res-file as list data structure'''
    def __init__(self, res_file):
        self.__resfilename = res_file
        self.__basefile = self.filename_wo_ending(self.__resfilename)


    def filename_wo_ending(self, resfilename):
        '''returns the input file name without ending'''
        file_ext = ''
        try:
            file_ext = os.path.splitext(resfilename)
            basefile = file_ext[0]
            if file_ext[1] != '.res':
                print("Please give a res-file name as argument!")
                sys.exit(0)
        except(AttributeError): 
            basefile = ''
        return str(basefile)

    
    def get_res_list(self):
        '''read the resfile and return it as list'''
        filename = self.__basefile+'.res'
        reslist = []
        try:
            with open(filename, 'r') as f:
                for line in f:
                    reslist.append(line)
        except(IOError):
            print('Unable to read {}'.format(filename))
            sys.exit()
        return reslist
        

    def write_resfile(self, reslist, ending):
        '''
        makes a copy of the resfile and writes basefile+ending.
        ending should be '.ins' or '.res'
        '''
#        try:
#            shutil.copyfile(self.__basefile+'.res', self.__basefile+'.bak')
#        except(IOError):
#            print 'Unable to create backup file!'
#            sys.exit(-1)
        try:
            #nfile = open(self.__basefile+'-new.res', 'w')  # open the res file and add '-new' to the name
            nfile = open(self.__basefile+ending, 'w')  # open the ins file
        except(IOError):
            print('Unable to write res file!')
            sys.exit(-1)
        for i in reslist:            #modified reslist
            nfile.write("%s" % i)    #write the new file
        nfile.close()
    
   # res_list = property(get_res_list, write_resfile)




class ResListEdit():
    '''Class to add, remove and list lines of the resfile
    '''
    def __init__(self, reslist, find_atoms):
        self.__reslist = reslist
        self._find_atoms = find_atoms

    
    def add_line(self, linenum, new_line):
        '''adds a single line of string to the res file
        checks are missing which make sure the data structure of new_line is correct!
        '''
        linenum = linenum-1
        # adds the new line:
        self.__reslist.insert(linenum, new_line+'\n')
        return self.__reslist

    
    def remove_line(self, linenum, rem=False, remove=False, frontspace=False):
        '''removes a single line from the res file with tree different methods. 
           The default is a space character in front of the line (frontspace). 
           This removes the line in the next refinement cycle. "rem" writes rem 
           in front of the line and "remove" clears the line.'''
        line = self.__reslist[linenum]
        if rem:   # comment out with 'rem ' in front
            self.__reslist[linenum] = 'rem '+line
            if misc.multiline_test(line):
                self.__reslist[linenum+1] = 'rem '+self.__reslist[linenum+1]
        elif remove:  # really delete the line "linenum"
            del self.__reslist[linenum]
            if misc.multiline_test(line):
                del self.__reslist[linenum]
        if frontspace:  # only put a space in front
            self.__reslist[linenum] = ' '+line
            if misc.multiline_test(line):
                self.__reslist[linenum+1] = ' '+self.__reslist[linenum+1]
        return self.__reslist

    
    
    def list_lines(self, startline, endline):
        '''returns the lines between startline and endline as list'''
        lines = []
        start = int(startline)
        end = int(endline)
        for line, i in enumerate(self.__reslist):
            if line >= startline:
                lines.append(i)
            if line == endline:
                break
        return lines
        
    
    def getAll(self):
        '''returns the whole resfile as list'''
        return self.__reslist

    
    def find_fvarlines(self):
        '''
        Finds the FVAR line or the first line with an atom.
        '''
        fvarlines = misc.find_multi_lines(self.__reslist, r'^FVAR.+[0-9]+')
        if not fvarlines:   # There is no FVAR in the res file after SHELXS!
            for num, i in enumerate(self.__reslist):
                if self._find_atoms.get_atom(i):
                    first_atom = num
                    break
            fvarlines = []
            fvarlines.append(first_atom-1)
            self.__reslist.insert(first_atom-1, ' \n')
        return fvarlines
    
    
    def insert_frag_fend_entry(self, atoms, fragline, fvarlines):
        '''
        Inserts the FRAG ... FEND entry in the res file.
        '''
        dblines = []
        db = atoms[:]
        db = [list(map(str, i)) for i in db]
        for i in db:
            i[2] = i[2].ljust(8, '0').rjust(9, ' ') 
            i[3] = i[3].ljust(8, '0').rjust(9, ' ')
            i[4] = i[4].ljust(8, '0').rjust(9, ' ')
            dblines.append('    '.join(i).rstrip())
        dblines = '\n'.join(dblines)        
        dblines = '  '.join(fragline)+'\n'+dblines
        dblines = '\n The following is from DSR:\n'+dblines
        dblines = dblines+'\nFEND\n\n'
        #dblines = misc.ll_to_string(db)+'\n'
        self.__reslist.insert(fvarlines[0]+2, dblines)   # insert the db entry right after FVAR


    def set_fvar(self, occupancynumber, fvarlines):
        '''
        Inserts additional free variables according to the occ parameter
        This function starts at the end of parse_dsr_line() so we don't have to care about it anywhere else. 
        '''
        occupancynumber = occupancynumber.strip('-')
        fvar_list = []
        
        for line in fvarlines:
            fvar = self.__reslist[line].split()
            if fvar:
                del fvar[0]
            fvar_list.extend(fvar)
        if len(fvar_list) != 0:
            for line in fvarlines:
                self.__reslist[line] = ' \n' # removes the old FVAR

        varlen = len(fvar_list) # how many numbers do we have?
        num = occupancynumber.split('.')   # the occupancynumber is split in the fvar part and the occupancy part
        fvar = int(num[0])//10              # e.g. 20.5 is fvar 2 and occupancy 0.5
        difference = (fvar - varlen)
        if difference > 0:
            for i in range(difference):
                fvar_list.append('0.5')   # if an fvar is missing, add a new one
        line_length = 7
        lines = []
                
        for i in range(0, len(fvar_list), line_length):
            l = 'FVAR  '+'  '.join( "{:<8}".format(x[:6].ljust(6, '0')) for x in fvar_list[i:i+line_length] )
            lines.append(l)
            fvars = '\n'.join(lines)
            fvars = fvars+'\n'
        self.__reslist.insert(fvarlines[0]+1, fvars)

    
    
# for testing
if __name__ == '__main__':
    from options import OptionsParser
    from dbfile import global_DB
    from atomhandling import FindAtoms
    from resfile import ResList, ResListEdit
    options = OptionsParser()
    res_list = ResList(options.res_file)
    reslist =  res_list.get_res_list()
    find_atoms = FindAtoms(reslist)
    rle = ResListEdit(reslist, find_atoms)
    fvarlines = rle.find_fvarlines()
    gdb = global_DB()
    db = gdb.build_db_dict()
    fragment = 'toluene'
    dbatoms = gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
    dbhead = gdb.get_head_from_fragment(fragment)        # this is only executed once
    
    fragline = gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
    rle.insert_frag_fend_entry(dbatoms, fragline, fvarlines)    
    rle.set_fvar('91.1234', fvarlines)
    

        
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
    

