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
import re
import sys
import misc
from resfile import ResListEdit
import constants
import textwrap


__metaclass__ = type  # use new-style classes

class DSR_Parser():
    '''
    handles the parsing of the DSR command
    
    This Class should have the reslist object as input and 
    should output the dictionary of the parsed dsr command.
    
    Additionally the FVAR line in the reslist gets corrected by set_fvar() at 
    the end of parse_dsr_line().
    '''
    def __init__(self, reslist, rle):
        self.__reslist = reslist
        self.__rle = rle
        self.__regex = '^rem\s{1,5}DSR\s{1,5}.*'.lower()
        self.__dsr_string = self.find_dsr_command(line=True).lower()
        self.__dsr = misc.makelist(self.__dsr_string)
    
    
    def find_dsr_command(self, line=False):
        '''
        line = False  -> Line number
        line = True  -> Text string
        find the lines with a DSR command entry and return its line number as default
        or the text string when line is set to True'''
        dsr_str = ''
        indexnum = [i for i, l in enumerate(self.__reslist) for m in [re.search(self.__regex, l.lower())] if m]
        try:
            indexnum[0]
        except(IndexError):
            print(' no proper DSR command found! \n\n Have you really saved your .res file?\n')
            sys.exit()

        if len(indexnum) > 1:
            print('only one dsr command at once is allowed! Exiting...')
            sys.exit(-1)
    
        if line is True:  # returns the string
            dsr_str = str(self.__reslist[indexnum[0]])
            if misc.multiline_test(dsr_str):
                # in case of a multiline command, strip the '=' and the newline
                dsr_str = dsr_str.rstrip('\n\r= ')
                dsr_str = re.sub(' +',' ', dsr_str)
                dsr_str = dsr_str+self.__reslist[indexnum[0]+1]
                return dsr_str                    # return the dsr text string
            else:
                # just return the dsr command as string in one line
                return dsr_str
        else:  # returns the index number of the command
            dsr_str = str(self.__reslist[indexnum[0]])
            if misc.multiline_test(dsr_str):
                # in case of a multiline command, strip the '=' and the newline
                dsr_str = dsr_str.rstrip('\n\r= ')
                dsr_str = dsr_str+' '+self.__reslist[indexnum[0]+1]
                dsr_list = dsr_str.split()
                txt = ' '.join(dsr_list)
                #for i in dsr_list:
                #    if i in constants.SHX_CARDS:
                #        print 'Error in DSR command line: Wrapped line with empty next line.'
                #        sys.exit()
                # wrap the line after 75 chars:
                dsrlines = '\n'.join(textwrap.wrap(txt, 75, initial_indent='rem ', subsequent_indent = 'rem ')) 
                dsrlines = dsrlines+'\n\n'
                self.__reslist[indexnum[0]] = '' # delete old line
                self.__reslist[indexnum[0]+1] = '' # delete old line
                self.__reslist.insert(indexnum[0]+1, dsrlines)
            else:
                dsr_str = re.sub(' +',' ', dsr_str)
                dsr_list = dsr_str.split()
               # dsr_list.insert(0, 'rem') # comment out the old line
                txt = ' '.join(dsr_list)
                dsrlines = '\n'.join(textwrap.wrap(txt, 75, initial_indent='rem ', subsequent_indent = 'rem ')) # wrap the line after 75 chars
                dsrlines = dsrlines+'\n\n'
                self.__reslist[indexnum[0]] = '' # delete old line
                self.__reslist.insert(indexnum[0]+1, dsrlines)
            if misc.multiline_test(dsr_str):
                return indexnum[0]+1
            else:
                return indexnum[0]                     # return the line index number
    

    def find_commands(self, command):
        '''returns the value of the input string argument as string'''
        # hier vielleicht sogar mit match, damit OCCC nicht gültig ist
        if command in self.__dsr:
            try:
                c_index = self.__dsr.index(command)+1  # find the index of the command+1
                cnum = self.__dsr[c_index]             # get the index value
            except(IndexError):
                cnum = False
        else:
            cnum = False
        return cnum


    def find_atoms(self, start, *stop):
        '''returns the source and target atoms between a single start argument
        and one or several stop arguments'''
        try:
            atindex = self.__dsr.index(start)+1
        except(ValueError):
            if start == 'WITH':
                print('No source atoms given!')
                sys.exit(-1)
            if start == 'ON':
                print('No target atoms given!')
                sys.exit(-1)
        atoms = []
        for i in self.__dsr[atindex:]: # start at the position of the first atom
            atoms.append(i.upper())
            if i in stop:   # stop if stop conditions met
                atoms.pop() # erase the stop condition from the list
                break
        return atoms


    def minimal_requirements(self):
        '''
        Checks if minimal requirements of the dsr command are met. E.g. the WITH and ON command 
        '''
        check = ('WITH', 'ON')
        for i in check:
            try:
                if i not in self.__dsr:
                    raise Exception
            except:
                print('\nNo "WITH" or "ON" statement in the dsr command line found!')
                sys.exit()


    def parse_dsr_line(self):
        '''returns the different parameters from the dsr command as dict
        It needs find_commands() and find_atoms() to parse the line.
        '''
        self.minimal_requirements()
        # syntax:
        # rem dsr put|add|replace fragment with source on target part xx AFIX 17x occ occupancy
        # the command:
        command = self.__dsr[2]
        # get the fragment:
        fragment = self.__dsr[3]
        # make sure the command is correct:
        command_list = ('PUT', 'REPLACE', 'ADD')
        if command not in command_list:
            print('No proper command string found in DSR command line!\n')#, self.__dsr
            sys.exit(-1)
        # Source and target atoms:
        # In paerenteses are one start und one to multiple stop conditions:
        source = self.find_atoms('WITH', 'ON')
        target = self.find_atoms('ON', 'PART', 'OCC', 'RESI', 'DFIX', '')
        if 'RESI' in self.__dsr:
            residue = self.find_atoms('RESI', 'PART', 'OCC', 'RESI', 'DFIX', '')
            # RESI is True but no residue returns -> only RESI in command line
            # hence, return that we at least want residues from the db
            if not residue:
                residue = 'DB'
        else:
            residue = ''
        dsr_dict = {
            'command': str(command),
            'fragment': str(fragment),
            'source': source, 
            'target': target,
            'part': self.find_commands('PART'),
            'dfix': self.find_commands('DFIX'),
            'occupancy': self.find_commands('OCC'), 
            'resi': residue
            };
        return dsr_dict
    
    
    @property
    def fragment(self):
        return self.parse_dsr_line()['fragment']
    
    @property
    def occupancy(self):
        return self.parse_dsr_line()['occupancy']
    
    @property
    def command(self):
        self.parse_dsr_line()['command']
        
    @property
    def source(self):
        self.parse_dsr_line()['source']
    
    @property
    def target(self):
        self.parse_dsr_line()['target']    
        
    @property
    def resi(self):
        self.parse_dsr_line()['resi']
    
    @property
    def dfix(self):
        self.parse_dsr_line()['dfix']
    
    @property
    def source(self):
        self.parse_dsr_line()['source']    
        
    
    





#for testing:
if __name__ == '__main__':
    from resfile import ResList, ResListEdit
    from options import OptionsParser 
    options = OptionsParser()
    rl = ResList(options.res_file)
    reslist = rl.get_res_list()
    rle = ResListEdit(reslist, options.res_file)
    #dsr_line = dsrp.parse_dsr_line()
    dsrp = DSR_Parser(reslist, rle)    
    #dsr_line = dsrp.parse_dsr_line()
    #print dsr_line
 #   dsr_str = dsrp.find_dsr_command(line=False)
 
 #   dsr_num = dsrp.find_dsr_command(line=False)
 #   
  
    
    print('\ncheck for FVAR:')
    for i in reslist:
        if i.startswith('FVAR'):
            print(i)
    print()
    
#    for i in dsr_line:
#        print i.ljust(9), '=', str(dsr_line.get(i)).ljust(11)
#    print 
   

    print('String:', dsrp.find_dsr_command(line=True))
    print('Zeile mit dbentry:\n', dsrp.find_dsr_command(line=False), '\n')
    

    
