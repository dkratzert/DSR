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
import textwrap


__metaclass__ = type  # use new-style classes


class DSRParser():
    """
    handles the parsing of the DSR command

    This Class should have the reslist object as input and
    should output the dictionary of the parsed dsr command.

    Additionally the FVAR line in the reslist gets corrected by set_fvar() at
    the end of parse_dsr_line().
    """
    def __init__(self, reslist):
        """
        :param reslist: list of strings of .res file
        """
        self._reslist = reslist
        self._dsr_regex = '^rem\s{1,5}DSR\s{1,5}.*'
        self._dsr_list = self.find_dsr_command(line=True)
        # Indicates if RESI is on of off in dsr command. Turns True if activated:
        # class and number are defined in self.resi
        self._resiflag = False
        try:
            self.dsr_dict = self.parse_dsr_line()
        except Exception:
            print("*** No valid DSR command line found. ***")
            raise

    def find_dsr_command(self, line=False):
        """
        line = False  -> Line number
        line = True  -> Text string
        find the lines with a DSR command entry and return its line number as
        default or the text string when line is set to True
        :param line: bool
        """
        hklf_endline = misc.find_line(self._reslist, r'^HKLF\s+[1-6]')
        dsr_str = ''
        multiline = False
        indexnum = misc.find_multi_lines(self._reslist, self._dsr_regex)
        try:
            line_number = int(indexnum[0])
        except IndexError:
            print('*** no proper DSR command found! \n'
                  'Have you really saved your .res file? ***\n')
            sys.exit()
        if int(line_number) > int(hklf_endline):
            print('*** A DSR command after HKLF is not allowed! '
                  'Check line {} ***'.format(line_number+1))
            sys.exit()
        if len(indexnum) > 1:
            print('*** Only one DSR command at once is allowed! ***')
            sys.exit(-1)
        # TODO: rip out dsr command removal:
        if line:  # returns the string
            dsr_str = str(self._reslist[line_number])
            if misc.multiline_test(dsr_str):
                multiline = True
                # in case of a multiline command, strip the '=' and the newline
                dsr_str = dsr_str.rstrip('\n\r= ')
                dsr_str = re.sub(' +', ' ', dsr_str)
                dsr_str = dsr_str+self._reslist[line_number+1]
                return dsr_str.upper().split()           # return the dsr command as list
            else:
                # just return the dsr command as string in one line
                return dsr_str.upper().split()
        else:  # returns the index number of the command
            dsr_str = str(self._reslist[line_number])
            dsr_str = re.sub(' +', ' ', dsr_str)
            if misc.multiline_test(dsr_str):
                multiline = True
                # in case of a multiline command, strip the '=' and the newline
                dsr_str = dsr_str.rstrip('\n\r= ')
                dsr_str = "{} {}".format(dsr_str, self._reslist[line_number+1])
            dsr_list = dsr_str.upper().split()
            txt = ' '.join(dsr_list)
            # wrap the line after 75 chars:
            dsrlines = textwrap.wrap(txt, 75, initial_indent='REM ', subsequent_indent='REM !')
            dsrlines = '\n'.join(dsrlines)
            dsrlines += '\n'
            self._reslist[line_number] = ''  # delete old line
            if multiline:
                self._reslist[line_number+1] = ''  # delete old line
            self._reslist[line_number] = dsrlines
            return line_number                     # return the line index number

    def find_commands(self, command):
        """
        returns the value of the input string argument as string
        :param command: string
        """
        if command in self._dsr_list:
            try:
                command_index = self._dsr_list.index(command)+1  # find the index of the command+1
                command_argument = self._dsr_list[command_index]             # get the index value
            except IndexError:
                command_argument = ''
        else:
            command_argument = ''
        return command_argument

    def find_atoms(self, start, *stop):
        """
        Returns the parameter between a single start argument
        and one or several stop arguments as list of strings. This can be atoms or
        other Parameters

        ----------
        :type start: str 
        :type stop: str
        :rtype: list
        """
        atindex = 0
        try:
            atindex = self._dsr_list.index(start)+1
        except ValueError:
            if start == 'WITH':
                print('*** No source atoms given! ***')
                sys.exit()
            if start == 'ON':
                print('*** No target atoms given! ***')
                sys.exit()
        atoms = []
        for i in self._dsr_list[atindex:]: # start at the position of the first atom
            atoms.append(i.upper())
            if i in stop:   # stop if stop conditions met
                atoms.pop() # erase the stop condition from the list
                break
        return atoms

    def minimal_requirements(self):
        """
        Checks if minimal requirements of the dsr command are met.
        E.g. the WITH and ON command
        """
        check = ('WITH', 'ON')
        for i in check:
            try:
                if i not in self._dsr_list:
                    raise Exception
            except:
                print('\n*** No "WITH" or "ON" statement in the dsr command line found! ***')
                sys.exit()

    def parse_dsr_line(self):
        """
        returns the different parameters from the dsr command as dict
        It needs find_commands() and find_atoms() to parse the line.

        >>> from resfile import ResList, ResListEdit
        >>> res_file = 'p21c.res'
        >>> rl = ResList(res_file)
        >>> reslist = rl.get_res_list()
        >>> rle = ResListEdit(reslist, res_file)
        >>> #dsr_line = dsrp.get_dsr_dict
        >>> dsrp = DSRParser(reslist)
        >>> dic = dsrp.all
        >>> l = sorted(dic)
        >>> for i in l:
        ...     print('{}: '.format(i), dic[i])
        command:  PUT
        dfix:  False
        fragment:  OC(CF3)3
        occupancy:  -31
        part:  2
        resi:  ['CF3']
        source:  ['O1', 'C1', 'C2', 'C3', 'C4']
        split:  False
        target:  ['O1_3', 'C1_3', 'Q6', 'Q4', 'Q7']
        """
        source = None
        if self.cf3_active:
            splitatom = False
            if 'SPLIT' in self._dsr_list:
                splitatom = True
        else:
            splitatom = False
            self.minimal_requirements()
        # syntax:
        # rem dsr put|add|replace fragment with source on target part xx AFIX 17x occ occupancy
        command = self._dsr_list[2]
        # get the fragment:
        fragment = self._dsr_list[3]
        if fragment in ['CF3', 'CF9'] and 'SPLIT' in self._dsr_list:
            print('*** Illegal combination of CF3 or CF9 with SPLIT! \nOnly CF6 with SPLIT is alowed! ***')
            sys.exit()
        # make sure the command is correct:
        command_list = ('PUT', 'REPLACE', 'ADD')
        if command not in command_list:
            print('*** No proper put/replace string found in DSR command line! ***')#, self._dsr_list
            sys.exit()
        for com in command_list:
            if com in self._dsr_list[3:]:
                print('*** PUT or REPLACE must occure once and only directly after REM DSR ***')
                sys.exit()
        # Source and target atoms:
        # In parenteses are one start und one to multiple stop conditions:
        if not self.cf3_active:
            # we need no soure atoms for cf3 groups
            source = self.find_atoms('WITH', 'ON')
        target = self.find_atoms('ON', 'PART', 'OCC', 'RESI', 'DFIX', 'SPLIT', '')
        if 'RESI' in self._dsr_list:
            residue = self.find_atoms('RESI', 'PART', 'OCC', 'RESI', 'DFIX', '')
            # RESI is in dsr_list but no residue returns -> only RESI in command line:
            if not residue:
                residue = []  # No residue defined, but RESI is active
            self._resiflag = True  # means residue is active
        else:
            residue = False
            self._resiflag = False  # residue completely turned off
        dfix = False
        if 'DFIX' in self._dsr_list:
            dfix = True
        part = self.find_commands('PART')
        if part:
            try:
                float(part)
            except ValueError:
                print('*** Part without numerical value supplied.',
                      'Please give a part number after PART in the DSR command. ***')
                sys.exit(False)
            if float(part) > 999:
                print('*** only 999 parts allowed in SHELXL! ***')
                sys.exit(False)
            try:
                if len(part) > 4:
                    print('*** Illegal part number supplied: {}. Please give just one digit '
                          '(positive or negative) in the DSR command ***'.format(part))
                    sys.exit(False)
            except TypeError:
                # no part specified
                pass
        occupancy = self.find_commands('OCC')
        badocc_message = '*** Occupancy without numerical value supplied. Please define occupancy value after OCC ***'
        badocc_status = False
        if occupancy:
            num = occupancy.split('.')
            fvar = abs(int(num[0]))//10
            try:
                if float(fvar) > 99:
                    print('*** Only 99 free variables allowed in SHELXL! ***')
                    sys.exit()
            except ValueError:
                badocc_status = True
        if 'OCC' in self._dsr_list and not occupancy:
            badocc_status = True
        if badocc_status:
            print(badocc_message)
            sys.exit()
        dsr_dict = {
            'command': str(command),
            'fragment': str(fragment),
            'source': source,
            'target': target,
            'part': part,
            'dfix': dfix,
            'occupancy': occupancy,
            'resi': residue,
            'split': splitatom
            }
        return dsr_dict

    @property
    def cf3_active(self):
        """
        Does it have the CF3-Group keyword?
        :return: True/False
        """
        cf = False
        if self._dsr_list[3] in ['CF3', 'CF6', 'CF9']:
            cf = True
        return cf

    @property
    def fragment(self):
        """
        database fragment name
        """
        return self.dsr_dict['fragment'].lower()

    @property
    def occupancy(self):
        """
        occupancy of the fragment
        """
        return self.dsr_dict['occupancy']

    @occupancy.setter
    def occupancy(self, value):
        """
        occupancy of the fragment
        """
        self.dsr_dict['occupancy'] = value

    @property
    def command(self):
        return self.dsr_dict['command']

    @property
    def part(self):
        return self.dsr_dict['part']

    @property
    def source(self):
        return self.dsr_dict['source']

    @property
    def target(self):
        return self.dsr_dict['target']
    
    @property
    def split(self):
        return self.dsr_dict['split']

    @split.setter
    def split(self, value):
        self.dsr_dict['split'] = value

    @property
    def resi(self):
        """
        resi: empty list, class, number or class and number
        :rtype: list
        """
        return self.dsr_dict['resi']

    @property
    def resiflag(self):
        """
        defines if residue is on or off
        :rtype: bool
        """
        return self._resiflag

    @resiflag.setter
    def resiflag(self, value):
        """ 
        Setter for resiflag
        :type value: bool
        """
        self._resiflag = value

    @property
    def dfix(self):
        """
        dfix: bool True/False
        """
        return self.dsr_dict['dfix']

    @property
    def all(self):
        """
        Returns the complete dsr command dictionary
        :return: dict 
        """
        return self.dsr_dict


if __name__ == '__main__':
    
    import doctest

    failed, attempted = doctest.testmod()  # verbose=True)
    if failed == 0:
        print('passed all {} tests!'.format(attempted))
    else:
        print('{} of {} tests failed'.format(failed, attempted))