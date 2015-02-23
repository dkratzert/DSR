#-*- encoding: utf-8 -*-
#mÃ¶p
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
from atomhandling import Elem_2_Sfac, rename_dbhead_atoms, FindAtoms,\
    get_atomtypes, NumberScheme, SfacTable
import misc
import os
import sys
import constants
import re
import fnmatch
from resfile import ResList, ResListEdit
from dsrparse import DSR_Parser
from dbfile import global_DB
from misc import ll_to_string

__metaclass__ = type  # use new-style classes


def write_dbhead_to_file(filename, dbhead, resi_class, resi_number):
    '''
    write the restraints to an external file
    :param filename:     filename of database file
    :param dbhead:       database header
    :param resi_class:   SHELXL residue class
    :param resi_number:  SHELXL residue number
    :return filename:    full file name where restraints will be written
    '''
    number = '1'
    files = []
    # find a unique number for the restraint file:
    for file in misc.sortedlistdir('.'):
        if fnmatch.fnmatch(file, 'dsr_*_'+filename):
            filenum = file.split('_')
            if str.isdigit(filenum[1]):
                files.append(filenum[1])
    try:
        number = str(int(files[-1])+1)
    except(IndexError):
        pass
    if not resi_number and not resi_class:   # no residues
        filename = 'dsr_'+number+'_'+filename
    if resi_number:                          # only residue number known
        filename = 'dsr_'+resi_class+'_'+resi_number+'_'+filename
    if not resi_number and resi_class:       # only residue class known
        filename = 'dsr_'+resi_class+'_'+filename
    if os.path.isfile(filename):
        print('Previous restraint file found.'\
            ' Using restraints from "{}"'.format(filename))
        return filename
    else:
        print('Restraints were written to "{}"'.format(filename))
    try:
        dfix_file = open(filename, 'w')  # open the ins file
    except(IOError):
        print('Unable to write res file! Check directory write permissions.')
        sys.exit(False)
    for i in dbhead:            #modified reslist
        dfix_file.write("%s" %i)    #write the new file
    dfix_file.close()
    return filename


class InsertAfix(object):
    '''
    methods for the AFIX entry
    - dbhead is modified by Resi() if residues are used!
      RESI num class ist inserted there
    '''

    def __init__(self, reslist, dbatoms, fragment_atom_types, dbhead, dsr_line_dict, sfac_table,
                find_atoms, numberscheme, dfix_head=False):
        '''
        :param reslist:      list of the .res file
        :type reslist: list
        :param dbatoms:      list of the atoms in the database entry
                             [['O1', 3, '-0.01453', '1.66590', '0.10966'], ... ]
        :type dbatoms: list
        :param fragment_atom_types:      ['N', 'C', 'C', 'C']
        :param dbhead:       database header
        :param dsr_line_dict: {'occupancy': 'str',
                                    'dfix': True/False,
                                    'part': 'str',
                                    'fragment': 'str',
                                    'resi': ['class'] or ['number'] or ['class', [number]],
                                    'target': ['atom1', 'atom2', ...],
                                    'source': ['atom1', 'atom2', ...],
                                    'command': 'PUT/REPLACE'}
        :param sfac_table:   SHELXL SFAC table as list like: ['C', 'H', 'O', 'F', 'Al', 'Ga']
        :param find_atoms:   FindAtoms() object
        :param numberscheme: atoms numbering scheme like: ['O1', 'C1', 'C2', 'F1', 'F2', 'F3', 'C3']^^
        '''
        self._reslist = reslist
        self._find_atoms = find_atoms
        self._dbatoms = dbatoms
        self._dbhead = dbhead
        self.dfix_head = dfix_head
        self._fragment_atom_types = fragment_atom_types
        self._sfac_table = sfac_table
        self.numberscheme = numberscheme
        self.part = dsr_line_dict['part']
        self.occ = dsr_line_dict['occupancy']
        self.source_atoms = dsr_line_dict['source']
        self.target_atoms = dsr_line_dict['target']
        self._dfix = dsr_line_dict['dfix']


    def insert_dsr_warning(self):
        '''
        information to insert into .res-file
        '''
        return 'rem the following was inserted by DSR:\n'


    def remove_duplicate_restraints(self, dbhead, residue_class=''):
        '''
        removes restraints from the header which are already
        in the res-file

        :param dbhead:         database header (list of strings)
        :param residue_class:  SHELXL residue class
        :type residue_class:   string
        '''
        modified = False
        newhead = dbhead[:]
        for resline in self._reslist:
            resline = resline.strip().split()
            for num, headline in enumerate(dbhead):
                headline = headline.strip().split()
                if not headline:
                    continue
                if headline == resline and headline[0][:4] in constants.RESTRAINT_CARDS:
                    # remove the restraint:
                    newhead[num] = '' #'rem '+newhead[num]
                    modified = True
                    break
        if modified:
            print('\nAlready existing restraints for residue "{}" were not '
                    'applied again.'.format(residue_class))
        return newhead


    def distance_and_other_restraints(self, dbhead):
        '''
        Devides header in distance restraints (distance)
        and all other lines (others)
        Restraints are instead inserted after fragment fit

        :param dbhead:  database header
        '''
        distance = []
        others = []
        for num, headline in enumerate(dbhead):  # @UnusedVariable
            headline = headline.strip().split()
            try:
                headline[0]
            except(IndexError):
                continue
            if headline[0][:4] in constants.DIST_RESTRAINT_CARDS:
                distance.append(' '.join(headline)+'\n')
            else:
                others.append(' '.join(headline)+'\n')
        return [distance, others]



    def build_afix_entry(self, external_restraints, dfx_file_name, resi):
        '''
        build an afix entry with atom coordinates from the target atoms

        :param external_restraints:  True/False decision if restraints should be
                                     written to external file
        :param dfx_file_name:        name of file for external restraints
        :param resi:                Residue() object
        '''
        afix_list = []   # the final list with atoms, sfac and coordinates
        e2s = Elem_2_Sfac(self._sfac_table)
        new_atomnames = list(reversed(self.numberscheme)) # i reverse it to pop() later
        if resi.get_residue_class:
            self._dbhead = self.remove_duplicate_restraints(self._dbhead, 
                                                            resi.get_residue_class)
            if not external_restraints:
                self._dbhead = self._dbhead+['RESI {} {}'.format(
                                        resi.get_residue_class, resi.get_resinumber)]
            self._dbhead = resi.format_restraints(self._dbhead)
            #print(self._dbhead)
        else:
            # applies new naming scheme to head:
            old_atoms = [ i[0] for i in self._dbatoms]
            self._dbhead = rename_dbhead_atoms(new_atomnames, old_atoms, self._dbhead)
        # decide if restraints to external file or internal:
        distance_and_other = self.distance_and_other_restraints(self._dbhead)
        distance = distance_and_other[0]
        other_head = distance_and_other[1]
        if external_restraints:
            # in case of dfix, write restraints to file after fragment fit
            self._dbhead = misc.wrap_headlines(distance)
            # returns the real name of the restraints file:
            if self.dfix_head:
                dfx_file_name = write_dbhead_to_file(dfx_file_name, self.dfix_head, 
                                    resi.get_residue_class, resi.get_resinumber)
            else:
                dfx_file_name = write_dbhead_to_file(dfx_file_name, self._dbhead, 
                                    resi.get_residue_class, resi.get_resinumber)
            if self.dfix_head:
                self._dbhead = other_head
            self._dbhead = self._dbhead+['RESI {} {}'.format(
                                        resi.get_residue_class, resi.get_resinumber)]
        else:
            if self.dfix_head:
                self._dbhead = other_head+self.dfix_head
            self._dbhead = misc.wrap_headlines(self._dbhead)
        # list of atom types in reverse order
        reversed_fragm_atom_types = list(reversed(self._fragment_atom_types))
        coordinates = self._find_atoms.get_atomcoordinates(self.target_atoms)
        #target = self.target_atoms[:]  # a copy because we edit it later
        # a list of zeroed atom coordinates (afix_list) is built:
        for i in self._dbatoms:
            l = []
            sfac_num = str(e2s.elem_2_sfac(reversed_fragm_atom_types.pop()))
            l.insert(0, str(i[0]))         # Atomname
            l.insert(1, sfac_num)  # SFAC number
            l.insert(2, '0       ')
            l.insert(3, '0       ')
            l.insert(4, '0       ')
            l.insert(5, '11.000')
            l.insert(6, '0.04')
            afix_list.append(l)
        # for every atoms both in afix list and source_atoms, change the
        # coordinates to the respective value of the target atom:
        ind = 0
        for n, i in enumerate(afix_list):
            if i[0].upper() in self.source_atoms:
                ind = self.source_atoms.index(i[0].upper())
                try:
                    afix_list[n][2:5] =  coordinates[self.target_atoms[ind]]
                except(IndexError):
                    print('More source than target atoms present! Exiting...')
                    sys.exit(False)
        newlist = []
        for i in afix_list:
            i[0] = new_atomnames.pop()
            i[2] = '{:>10.6f}'.format(float(i[2]))
            i[3] = '{:>10.6f}'.format(float(i[3]))
            i[4] = '{:>10.6f}'.format(float(i[4]))
            newlist.append('   '.join(i).rstrip())
        atoms = '\n'.join(newlist)
        afixnumber = '179'   # makes afix 179 default
        if not self.occ:
            self.occ = ''
        if self.part:
            part = 'PART '+str(self.part)+' '+str(self.occ)+'\n'
            part2 = 'PART 0\n'
        else:
            part = ''
            part2 = ''
        if resi.get_residue_class:
            resinum = 'RESI 0\n'
        else:
            resinum = ''
        if external_restraints:
            if resi.get_residue_class:
                self._dbhead.append('\nREM The restraints for residue {} are in this'\
                    ' file:\n+{}\n'.format(resi.get_residue_class, dfx_file_name))
            else:
                self._dbhead.append('\nREM The restraints for this moiety are in this'\
                          ' file:\n+{}\n'.format(dfx_file_name))
        self._dbhead = ''.join(self._dbhead)
        warn = self.insert_dsr_warning()
        afix = warn+self._dbhead+part+'AFIX '+str(afixnumber)+'\n'+atoms+(
                '\nAFIX 0\n'+part2+resinum+'rem The end of the DSR entry\n\n')
        #print(self._dbhead)
        return afix


if __name__ == '__main__':
    # from resi import Resi
    res_file = 'p21c.res'
    invert = True
    rl = ResList(res_file)
    reslist = rl.get_res_list()
    dsrp = DSR_Parser(reslist, rl)
    dsr_dict = dsrp.parse_dsr_line()
    find_atoms = FindAtoms(reslist)
    rle = ResListEdit(reslist, find_atoms)
    gdb = global_DB(invert)
    db = gdb.build_db_dict()
    fragment = 'OC(CF3)3'
    fragline = gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
    dbatoms = gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
    dbhead = gdb.get_head_from_fragment(fragment)        # this is only executed once
    print(dbhead)
    resi = True #gdb.get_resi_from_fragment(fragment)
    dbtypes = get_atomtypes(dbatoms)
    #resi = Resi(reslist, dsr_dict, dbhead, residue, find_atoms)
    #dbhead = resi.make_resihead()

    sf = SfacTable(reslist, dbtypes)
    sfac_table = sf.set_sfac_table()
    num = NumberScheme(reslist, dbatoms, resi)
    numberscheme = num.get_fragment_number_scheme()


    afix = InsertAfix(reslist, dbatoms, dbtypes, dbhead, dsr_dict, sfac_table, find_atoms, numberscheme)
    print(afix.build_afix_entry(False, res_file, 'CF3'))




