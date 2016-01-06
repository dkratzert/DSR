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
import fnmatch
from resfile import ResList, ResListEdit
from dsrparse import DSR_Parser
from dbfile import global_DB
from constants import RESTRAINT_CARDS
from resi import Resi


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
    for filen in misc.sortedlistdir('.'):
        if fnmatch.fnmatch(filen, 'dsr_*_'+filename):
            filenum = filen.split('_')
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

def add_residue_to_dfix(dfix_head, resinum):
    '''
    Add a residue to a list of DFIX/DANG restraints
    DFIX 1.234 C1 C2 -> DFIX 1.234 C1_4 C2_4 
    >>> add_residue_to_dfix(['DFIX 1.456 C1 C2', 'DFIX 1.212 C3 C4'], 4)
    ['DFIX  1.456  C1_4  C2_4\\n', 'DFIX  1.212  C3_4  C4_4\\n']
    >>> add_residue_to_dfix(['DFIX 1.456 C1 C2', 'DFIX 1.212 C3 C4'], '5')
    ['DFIX  1.456  C1_5  C2_5\\n', 'DFIX  1.212  C3_5  C4_5\\n']
    '''
    newhead = []
    for line in dfix_head:
        line = line.split()
        try:
            line[3]
        except:
            newhead.append(' '.join(line))
            continue
        line[2] = line[2]+'_'+str(resinum)
        line[3] = line[3]+'_'+str(resinum)
        line = '  '.join(line)+'\n'
        newhead.append(line)
    return newhead

class InsertAfix(object):
    '''
    methods for the AFIX entry
    - dbhead is modified by Resi() if residues are used!
      RESI num class ist inserted there
    '''

    def __init__(self, reslist, dbatoms, fragment_atom_types, dbhead, dsr_line_dict, sfac_table,
                find_atoms, numberscheme, options, dfix_head=False):
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
        :param numberscheme: atoms numbering scheme like: ['O1', 'C1', 'C2', 'F1', 'F2', 'F3', 'C3']

        >>> res_file = 'p21c.res'
        >>> invert = True
        >>> rl = ResList(res_file)
        >>> reslist = rl.get_res_list()
        >>> dsrp = DSR_Parser(reslist, rl)
        >>> dsr_dict = dsrp.get_dsr_dict
        >>> find_atoms = FindAtoms(reslist)
        >>> rle = ResListEdit(reslist, find_atoms)
        >>> gdb = global_DB(invert)
        >>> db = gdb.build_db_dict()
        >>> fragment = 'pph3'
        >>> fragline = gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
        >>> dbatoms = gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
        >>> dbhead = gdb.get_head_from_fragment(fragment)        # this is only executed once
        >>> print(dbhead)
        ['SADI C1 C2 C2 C3 C3 C4 C4 C5 C5 C6 C6 C1 C7 C8 C8 C9 C9 C10 C10 C11 C11 C12 C12 C7 C13 C14 \
C14 C15 C15 C16 C16 C17 C17 C18 C18 C13', 'SADI P1 C1 P1 C7 P1 C13', \
'SADI 0.04 C1 C3 C1 C5 C5 C3 C4 C2 C4 C6 C2 C6 C7 C9 C7 C11 C9 C11 C10 C8 C10 C12 C8 C12 \
C13 C15 C13 C17 C15 C17 C16 C14 C16 C18 C14 C18', 'SADI 0.04 P1 C14 P1 C18 P1 C12 P1 C8 P1 \
C2 P1 C6', 'FLAT C1 > C6 P1', 'FLAT C7 > C12 P1', 'FLAT C13 > C18 P1', 'SIMU P1 > C18', \
'RIGU P1 > C18']
        >>> basefilename = filename_wo_ending(res_file)
        >>> resi = True #gdb.get_resi_from_fragment(fragment)
        >>> dbtypes = get_atomtypes(dbatoms)
        >>> #resi = Resi(reslist, dsr_dict, dbhead, residue, find_atoms)
        >>> resi = Resi(reslist, dsr_dict, dbhead='RESI CF3', db_residue_string='CF3', find_atoms=find_atoms)
        No residue number was given. Using residue number 4.
        >>> sf = SfacTable(reslist, dbtypes)
        >>> sfac_table = sf.set_sfac_table()
        >>> num = NumberScheme(reslist, dbatoms, resi)
        >>> numberscheme = num.get_fragment_number_scheme()
        RESI instruction is enabled. Leaving atom numbers as they are.
        >>> dfix_head = ''
        >>> afix = InsertAfix(reslist, dbatoms, dbtypes, dbhead, dsr_dict, sfac_table, find_atoms, numberscheme, options, dfix_head)
        >>> afix.remove_duplicate_restraints(dbhead, 'PPh3')
        ['SADI C1 C2 C2 C3 C3 C4 C4 C5 C5 C6 C6 C1 C7 C8 C8 C9 C9 C10 C10 C11 C11 C12 C12 C7 \
C13 C14 C14 C15 C15 C16 C16 C17 C17 C18 C18 C13', 'SADI P1 C1 P1 C7 P1 C13', \
'SADI 0.04 C1 C3 C1 C5 C5 C3 C4 C2 C4 C6 C2 C6 C7 C9 C7 C11 C9 C11 C10 C8 C10 C12 C8 C12 \
C13 C15 C13 C17 C15 C17 C16 C14 C16 C18 C14 C18', 'SADI 0.04 P1 C14 P1 C18 P1 C12 P1 C8 P1 \
C2 P1 C6', 'FLAT C1 > C6 P1', 'FLAT C7 > C12 P1', 'FLAT C13 > C18 P1', 'SIMU P1 > C18', \
'RIGU P1 > C18']
        >>> afix.build_afix_entry(False, basefilename+'.dfix', resi)
        'rem the following was inserted by DSR:\\nSADI_CF3 C1 C2 C2 C3 C3 C4 C4 C5 C5 C6 C6 \
C1 C7 C8 C8 C9 C9 C10 C10 C11 C11 =\\n   C12 C12 C7 C13 C14 C14 C15 C15 C16 C16 C17 C17 C18 \
C18 C13\\nSADI_CF3 P1 C1 P1 C7 P1 C13\\nSADI_CF3 0.04 C1 C3 C1 C5 C5 C3 C4 C2 C4 C6 C2 C6 C7 \
C9 C7 C11 C9 C11 C10 C8 =\\n   C10 C12 C8 C12 C13 C15 C13 C17 C15 C17 C16 C14 C16 C18 C14 \
C18\\nSADI_CF3 0.04 P1 C14 P1 C18 P1 C12 P1 C8 P1 C2 P1 C6\\nFLAT_CF3 C1 > C6 P1\\nFLAT_CF3 \
C7 > C12 P1\\nFLAT_CF3 C13 > C18 P1\\nSIMU_CF3 P1 > C18\\nRIGU_CF3 P1 > C18\\nRESI CF3 \
4\\nPART 2 -31\\nAFIX 179\\n\
P1   7     0.000000     0.000000     0.000000   11.000   0.04\\n\
C1   1     0.033372     0.232537     0.337248   11.000   0.04\\n\
C2   1     0.024600     0.306100     0.306300   11.000   0.04\\n\
C3   1     0.122500     0.190800     0.297100   11.000   0.04\\n\
C4   1    -0.109100     0.204000     0.333000   11.000   0.04\\n\
C5   1     0.000000     0.000000     0.000000   11.000   0.04\\n\
C6   1     0.000000     0.000000     0.000000   11.000   0.04\\n\
C7   1     0.000000     0.000000     0.000000   11.000   0.04\\n\
C8   1     0.000000     0.000000     0.000000   11.000   0.04\\n\
C9   1     0.000000     0.000000     0.000000   11.000   0.04\\n\
C10   1     0.000000     0.000000     0.000000   11.000   0.04\\n\
C11   1     0.000000     0.000000     0.000000   11.000   0.04\\n\
C12   1     0.000000     0.000000     0.000000   11.000   0.04\\n\
C13   1     0.000000     0.000000     0.000000   11.000   0.04\\n\
C14   1     0.000000     0.000000     0.000000   11.000   0.04\\n\
C15   1     0.000000     0.000000     0.000000   11.000   0.04\\n\
C16   1     0.000000     0.000000     0.000000   11.000   0.04\\n\
C17   1     0.000000     0.000000     0.000000   11.000   0.04\\n\
C18   1     0.000000     0.000000     0.000000   11.000   0.04\\n\
AFIX 0\\nPART 0\\nRESI 0\\nrem The end of the DSR entry\\n\\n'
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
        self.options = options


    def insert_dsr_warning(self):
        '''
        information to insert into .res-file
        '''
        return 'rem the following was inserted by DSR:\n'

    def collect_all_restraints(self):
        '''
        collects all restraints in the resfile and returns a list with them
        [['RIGU_CF3', 'O1', '>', 'F9'], '...']
        '''
        all_restraints =[]
        for n, resline in enumerate(self._reslist):
            resline = resline.strip(' \n\r')
            resline = resline.split()
            try:
                resline[0][:4]
            except:
                continue
            if resline[0][:4] in RESTRAINT_CARDS:
                # see for the next four lines if the lines continues with "=":
                line = 0
                while resline[-1] == '=':
                    resline = resline[:-1]+self._reslist[n+line+1].split()
                    line = line + 1
                    if not resline[-1] == '=':
                        break
                    if line > 500:
                        break
                all_restraints.append(resline)
        return all_restraints
    
    def remove_duplicate_restraints(self, dbhead, residue_class=''):
        '''
        removes restraints from the header which are already
        in the res-file.

        :param dbhead:         database header (list of strings)
        :param residue_class:  SHELXL residue class
        :type residue_class:   string
        '''
        all_restraints = self.collect_all_restraints()
        modified = False
        newhead = dbhead[:]
        for num, headline in enumerate(dbhead):
            headline = headline.split()
            for restr in all_restraints:
                if headline == restr:
                    newhead[num] = ''
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
            self._dbhead = resi.format_restraints(self._dbhead)
            self._dbhead = self.remove_duplicate_restraints(self._dbhead, 
                                                            resi.get_residue_class)
            if not external_restraints:
                self._dbhead = self._dbhead+['RESI {} {}'.format(
                                        resi.get_residue_class, resi.get_resinumber)]
        else:
            # applies new naming scheme to head:
            old_atoms = [ i[0] for i in self._dbatoms]
            self._dbhead = rename_dbhead_atoms(new_atomnames, old_atoms, self._dbhead)
        # decide if restraints to external file or internal:
        distance_and_other = self.distance_and_other_restraints(self._dbhead)
        distance = distance_and_other[0]
        other_head = distance_and_other[1]
        if external_restraints and not self.options.rigid_group:
            # in case of dfix, write restraints to file after fragment fit
            self._dbhead = misc.wrap_headlines(distance)
            # returns the real name of the restraints file:
            if self.dfix_head:
                if resi.get_residue_class:
                    self.dfix_head = add_residue_to_dfix(self.dfix_head, resi.get_resinumber)
                dfx_file_name = write_dbhead_to_file(dfx_file_name, self.dfix_head, 
                                    resi.get_residue_class, resi.get_resinumber)
            else:
                dfx_file_name = write_dbhead_to_file(dfx_file_name, self._dbhead, 
                                    resi.get_residue_class, resi.get_resinumber)
                self._dbhead = self._dbhead = other_head
            if self.dfix_head:
                self._dbhead = other_head
            self._dbhead = self._dbhead+['RESI {} {}'.format(
                                        resi.get_residue_class, resi.get_resinumber)]
        else:
            if self.dfix_head:
                if resi.get_residue_class:
                    self.dfix_head = add_residue_to_dfix(self.dfix_head, resi.get_resinumber)
                self._dbhead = other_head+self.dfix_head
                if not resi.get_residue_class:
                    self._dbhead = rename_dbhead_atoms(new_atomnames, old_atoms, self._dbhead)
            self._dbhead = misc.wrap_headlines(self._dbhead)
        # list of atom types in reverse order
        reversed_fragm_atom_types = list(reversed(self._fragment_atom_types))
        coordinates = self._find_atoms.get_atomcoordinates(self.target_atoms)
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
        if external_restraints and not self.options.rigid_group:
            if resi.get_residue_class:
                self._dbhead.append('\nREM The restraints for residue {} are in this'\
                    ' file:\n+{}\n'.format(resi.get_residue_class, dfx_file_name))
            else:
                self._dbhead.append('\nREM The restraints for this moiety are in this'\
                          ' file:\n+{}\n'.format(dfx_file_name))
        if self.options.rigid_group:
            if resi.get_residue_class:
                self._dbhead = ''
                self._dbhead = self._dbhead+'RESI {} {}\n'.format(
                                        resi.get_residue_class, resi.get_resinumber)
            else:
                self._dbhead = ''
        else:
            self._dbhead = ''.join(self._dbhead)
        #warn = self.insert_dsr_warning()
        afix = self._dbhead+part+'AFIX '+str(afixnumber)+'\n'+atoms+(
                '\nAFIX 0\n'+part2+resinum+'\n\n')#+'rem The end of the DSR entry\n\n')
        #print(self._dbhead)
        return afix


if __name__ == '__main__':
    import sys
    import doctest
    failed, attempted = doctest.testmod()#verbose=True)
    if failed == 0:
        print('passed all {} tests!'.format(attempted))
        
###################################################################        
    sys.exit()
    res_file = 'p21c.res'
    invert = True
    rl = ResList(res_file)
    reslist = rl.get_res_list()
    dsrp = DSR_Parser(reslist, rl)
    dsr_dict = dsrp.get_dsr_dict
    find_atoms = FindAtoms(reslist)
    rle = ResListEdit(reslist, find_atoms)
    gdb = global_DB(invert)
    db = gdb.build_db_dict()
    fragment = 'pph3'
    fragline = gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
    dbatoms = gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
    dbhead = gdb.get_head_from_fragment(fragment)        # this is only executed once
    print(dbhead)

    resi = True #gdb.get_resi_from_fragment(fragment)
    dbtypes = get_atomtypes(dbatoms)
    #resi = Resi(reslist, dsr_dict, dbhead, residue, find_atoms)
    #dbhead = resi.make_resihead()
    resi = Resi(reslist, dsr_dict, dbhead, 'RESI PPH3', find_atoms)
    sf = SfacTable(reslist, dbtypes)
    sfac_table = sf.set_sfac_table()
    num = NumberScheme(reslist, dbatoms, resi)
    numberscheme = num.get_fragment_number_scheme()





