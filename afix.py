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
    get_atomtypes, NumberScheme
import misc
import os
import sys
import constants

__metaclass__ = type  # use new-style classes


def write_dbhead_to_file(filename, dbhead, resi_class, resi_number):
    '''
    write the restraints to an external file
    :param filename:     filename of database file
    :param dbhead:       database header
    :param resi_class:   SHELXL residue class
    :param resi_number:  SHELXL residue number
    '''
    if resi_number:
        filename = resi_class+'_'+resi_number+'_'+filename
    else:
        filename = resi_class+'_'+filename
    if os.path.isfile(filename):
        print('Previous restraint file found.'\
            ' Using restraints from "{}"'.format(filename))
        return filename
    else:
        print('Restraints were written to "{}"'.format(filename))
    try:
        dfix_file = open(filename, 'w')  # open the ins file
    except(IOError):
        print('Unable to write res file!')
        sys.exit(-1)
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
                find_atoms, numberscheme):
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
        self._fragment_atom_types = fragment_atom_types
        self._sfac = sfac_table
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


    def remove_duplicate_restraints(self, dbhead, residue_class):
        '''
        removes restraints from the header which are already
        in the res-file

        :param dbhead:         database header (list of strings)
        :param residue_class:  SHELXL residue class
        '''
        modified = False
        newhead = dbhead[:]
        for resline in self._reslist:
            resline = resline.strip().split()
            for num, headline in enumerate(dbhead):
                headline = headline.strip().split()
                if headline == resline and headline[0][:4] in constants.RESTRAINT_CARDS:
                    # remove the restraint:
                    newhead[num] = '' #'rem '+newhead[num]
                    modified = True
                    break
        if modified:
            print('\nAlready existing restraints for residue "{}" were not '
                    'applied again.'.format(residue_class))
        return newhead


    def remove_all_restraints(self, dbhead):
        '''
        Devides header in distance restraints (distance)
        and all other lines (oldhead)
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



    def build_afix_entry(self, external_restraints, dfx_file_name, residue_class):
        '''
        build an afix entry with atom coordinates from the target atoms

        :param external_restraints:  True/False decision if restraints should be
                                     written to external file
        :param dfx_file_name:             name of file for external restraints
        :param residue_class:        SHELXL residue class
        '''
        afix_list = []   # the final list with atoms, sfac and coordinates
        e2s = Elem_2_Sfac(self._sfac)
        new_atomnames = list(reversed(self.numberscheme)) # i reverse it to pop() later
        # all non-atoms between start tag and FRAG card with new names:
        dbhead = self._dbhead
        if residue_class:
            dbhead = self.remove_duplicate_restraints(dbhead, residue_class)
        else:
            # applies new naming scheme
            old_atoms = [ i[0] for i in self._dbatoms]
            dbhead = rename_dbhead_atoms(new_atomnames, old_atoms, dbhead)
        removed_restr = self.remove_all_restraints(dbhead)
        dbhead_for_dfix = removed_restr[0]
        dbhead_others = misc.wrap_headlines(removed_restr[1])
        if self._dfix:
            dbhead = dbhead_others
        if external_restraints and not self._dfix:
            # in case of dfix, write restraints ti file after fragment fit
            dbhead = dbhead_others
            resinumber = False
            dbhead_for_dfix = misc.wrap_headlines(dbhead_for_dfix)
            # returns the real name of the restraints file:
            dfx_file_name = write_dbhead_to_file(dfx_file_name, dbhead_for_dfix, residue_class, resinumber)
        if not external_restraints and not self._dfix:
            dbhead_for_dfix = misc.wrap_headlines(dbhead_for_dfix)
            dbhead = dbhead_others+dbhead_for_dfix
        # list of atom types in reverse order
        reversed_fragm_atom_types = list(reversed(self._fragment_atom_types))
        coordinates = self._find_atoms.get_atomcoordinates(self.target_atoms)
        #target = self.target_atoms[:]  # a copy because we edit it later
        # a list of zeroed atom coordinates (afix_list) is built:
        for i in self._dbatoms:
            l = []
            l.insert(0, str(i[0]))         # Atomname
            l.insert(1, str(e2s.elem_2_sfac(reversed_fragm_atom_types.pop())))  # SFAC number
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
                afix_list[n][2:5] =  coordinates[self.target_atoms[ind]]
        newlist = []
        for i in afix_list:
            i[0] = new_atomnames.pop()
            i[2] = '{:>10.6f}'.format(float(i[2]))
            i[3] = '{:>10.6f}'.format(float(i[3]))
            i[4] = '{:>10.6f}'.format(float(i[4]))
            newlist.append('   '.join(i).rstrip())
        atoms = '\n'.join(newlist)
        self.afixnumber = '179'   # makes afix 179 default
        if not self.occ:
            self.occ = ''
        if self.part:
            part = 'PART '+str(self.part)+' '+str(self.occ)+'\n'
            part2 = 'PART 0\n'
        else:
            part = ''
            part2 = ''
        if residue_class:
            resi = 'RESI 0\n'
        else:
            resi = ''
        if external_restraints and not self._dfix:
            dbhead.append('REM The restraints for residue_class {} are in this'\
                          ' file:\n+{}\n'.format(residue_class, dfx_file_name))
        dbhead = ''.join(dbhead)
        warn = self.insert_dsr_warning()
        afix = warn+dbhead+str(part)+'AFIX '+str(self.afixnumber)+'\n'+atoms+(
                '\nAFIX 0\n'+part2+resi+'rem The end of the DSR entry\n\n')
        return afix


if __name__ == '__main__':
    from dbfile import global_DB
    from dsrparse import DSR_Parser
    from atomhandling import SfacTable
    from resfile import ResList, ResListEdit
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
    fragment = 'PFAnion'
    fragline = gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
    dbatoms = gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
    dbhead = gdb.get_head_from_fragment(fragment)        # this is only executed once
    resi = True #gdb.get_resi_from_fragment(fragment)
    dbtypes = get_atomtypes(dbatoms)
    #resi = Resi(reslist, dsr_dict, dbhead, residue, find_atoms)
    #dbhead = resi.make_resihead()

    sf = SfacTable(reslist, dbtypes, res_file)
    sfac_table = sf.set_sfac_table()
    num = NumberScheme(reslist, dbatoms, resi)
    numberscheme = num.get_fragment_number_scheme()


    afix = InsertAfix(reslist, dbatoms, dbtypes, dbhead, dsr_dict, sfac_table, find_atoms, numberscheme)
    print(afix.build_afix_entry(True, False, False))




