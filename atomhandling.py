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
import re, sys
import string
from atoms import Element
from misc import find_line, find_multi_lines,\
    atomic_distance
from constants import atomregex, SHX_CARDS
from atoms import atoms


__metaclass__ = type  # use new-style classes


def get_atoms(atlist):
    """
    returns all atoms found in the input as list of lists
    >>> get_atoms(['F8    4    0.349210   0.073474   0.519443 -21.00000   0.03106', 'foo'])
    [['F8', '4', '0.349210', '0.073474', '0.519443']]
    """
    atoms = []
    try:
        atlist[0]
    except:
        return []
    if isinstance(atlist[0], list):
        for i in atlist:
            atoms.append(i[:5])
        return atoms
    for i in atlist:
        if re.match(atomregex, str(i)):  # search atoms
            l = i.split()[:5]  # convert to list and use only first 5 columns
            if l[0].upper() not in SHX_CARDS:  # exclude all non-atom cards
                atoms.append(l)
    return atoms


def get_atomtypes(dbatoms):
    """
    find all atoms in a list of shelxl format atom lines.
    returns a list like ['N', 'C', 'C', 'C'].

    :param dbatoms:  atoms in a database entry
                     [['O1', 3, '-0.01453', '1.66590', '0.10966'],
                      ['C1', 1, '-0.00146', '0.26814', '0.06351'],
                      [...]
                     ]
    :type dbatoms: list
    >>> atoms = [['O1', 3, '-0.01453', '1.66590', '0.10966'],\
                 ['C1', 1, '-0.00146', '0.26814', '0.06351'],\
                 ['Au15', 8, '0.00146', '0.26814', '0.06351']]
    >>> get_atomtypes(atoms)
    ['O', 'C', 'AU']
    >>> atoms = [['O1', 3, '-0.01453', '1.66590', '0.10966'],\
                 ['C1', 1, '-0.00146', '0.26814', '0.06351'],\
                 ['Ax15', 8, '0.00146', '0.26814', '0.06351']]
    >>> get_atomtypes(atoms)
    Traceback (most recent call last):
    ...
    KeyError
    """
    el = Element()
    found = []
    # find lines with atoms and see if they are in the atom list
    for i in dbatoms:
        sfacnum = i[1]
        atom_name = i[0].upper()    # i is the full atom name with number suffix like C1
        try:
            if int(sfacnum) < 0:
                found.append(el.get_element(abs(int(sfacnum))))
                continue
        except:
            pass
        found.append(el.get_atomlabel(atom_name))
    if len(dbatoms) != len(found):    # do we really need this here??
        print("*** One of the Atoms in the database entry is not correct! ***")
        raise KeyError
    return found


def replace_after_fit(rl, reslist, resi, fragment_numberscheme, cell):
    """
    deletes the atoms in replace mode that are near the fragment atoms

    :param rl: Reslist() instance
    :param reslist: .res file list
    :param resi: Resi() instance
    :param fragment_numberscheme: atom names of the fitting fragment
    :param cell: cell parameters

    >>> from resfile import ResList
    >>> from dbfile import global_DB
    >>> from resfile import ResListEdit
    >>> from dsrparse import DSR_Parser
    >>> from resi import Resi
    >>> res_file = 'p21c.res'
    >>> rl = ResList(res_file)
    >>> reslist = rl.get_res_list()
    >>> find_atoms = FindAtoms(reslist)
    >>> rle = ResListEdit(reslist, find_atoms)
    >>> dsrp = DSR_Parser(reslist, rle)
    >>> dsr_dict = dsrp.get_dsr_dict
    >>> dbhead = ['RESI CF3']
    >>> db_residue_string='CF3'
    >>> resi = Resi(reslist, dsr_dict, dbhead, db_residue_string, find_atoms)
    No residue number was given. Using residue number 4.
    >>> fragment_numberscheme = ['O1_1', 'C1_1', 'C2_1', 'F1_1', 'F2', 'F3', 'C3', 'F4', 'F5', 'F6', 'C4', 'F7', 'F8', 'F9']
    >>> cell = rle.get_cell()

    #>>> replace_after_fit(rl, reslist, resi, fragment_numberscheme, cell)
    """
    remdist=1.3
    from resfile import ResListEdit
    find_atoms = FindAtoms(reslist)
    if resi.get_resinumber:
        frag_at = []
        for i in fragment_numberscheme:
            at = i + '_{}'.format(resi.get_resinumber)
            frag_at.append(at)
    else:
        frag_at = fragment_numberscheme
    atoms_to_delete = find_atoms.find_atoms_to_replace(frag_at, cell, remdist)
    if atoms_to_delete:
        print('Replacing following atoms (< {0} A near fragment):\n'.format(remdist),
              ' '.join(sorted(set(atoms_to_delete))))
    target_lines = set(find_atoms.get_atom_line_numbers(atoms_to_delete))
    rle = ResListEdit(reslist, find_atoms)
    for i in target_lines:
        i = int(i)
        rle.remove_line(i, rem=False, remove=True, frontspace=False)
    rl.write_resfile(reslist, '.res')
    reslist = rl.get_res_list()
    return reslist, find_atoms


class FindAtoms():
    """
    Finds Atoms and creates a data structure
    """

    def __init__(self, reslist):
        '''
        :param reslist: SHELXL .res file as list
        '''
        self._reslist = reslist
        sfacline = find_multi_lines(self._reslist, r'SFAC\s+[a-zA-Z]+')
        try:
            self.sfac = self._reslist[sfacline[0]].split()[1:]
        except IndexError:
            print("*** No SFAC card found in .res file ***")
            sys.exit()
        self.e2s = Elem_2_Sfac(self.sfac)
        self._residues = self.collect_residues()

    @property
    def atoms_as_residues(self):
        """
        returns   residues = {   #  0           1            2         3      4      5          6
                            {'0': ['C1', ['x', 'y', 'z'], linenumber, class, part, element, sfac_number],
                           ['C2', ['x', 'y', 'z'], linenumber, class, part, element, sfac_number]},
                     {'1': ['C1', ['x', 'y', 'z'], linenumber, class, part, element, sfac_number],
                     []} }
        """
        return self._residues

    @staticmethod
    def is_atom(atomline):
        """
        returns all atoms found in the input as list if they are real atoms

        :param atomline:  'O1    3    0.120080   0.336659   0.494426  11.00000   0.01445 ...'
        >>> FindAtoms.is_atom(atomline = 'O1    3    0.120080   0.336659   0.494426  11.00000   0.01445 ...')
        ['O1', '3', '0.120080', '0.336659', '0.494426']
        >>> FindAtoms.is_atom(atomline = 'O1    0.120080   0.336659   0.494426  11.00000   0.01445 ...')
        []
        >>> FindAtoms.is_atom(atomline = 'O1    0.120080   0.336659   0.494426  11.00000   0.01445 ...')
        []
        """
        if re.search(atomregex, str(atomline)):        # search atoms
            atom = atomline.split()[:5]              # convert to list and use only first 5 columns
            if atom[0].upper() not in SHX_CARDS:      # exclude all non-atom cards
                return atom
            else:
                return []
        else:
            return []

    @staticmethod
    def get_resi_definition_dict(resi):
        """
        returns the residue number and class of a string like 'RESI TOL 1'
        or 'RESI 1 TOL'
        {'class': 'TOL', 'number': '1'}

        :param resi: ['RESI', 'number', 'class']
        :type resi: list or string

        >>> sorted(list(FindAtoms.get_resi_definition_dict('RESI 1 TOL'))
        ['class', 'number']
        >>> sorted(FindAtoms.get_resi_definition_dict('RESI 1 TOL').values())
        ['1', 'TOL']
        >>> FindAtoms.get_resi_definition_dict('RESI 1 TOL')
        {'class': 'TOL', 'number': '1'}
        """
        resi_dict = {
            'class' : None,
            'number': None}
        try:
            resi.remove('RESI')
        except AttributeError:
            resi = resi.split()
            resi.remove('RESI')
        resi.sort()
        if str.isalpha(resi[-1][0]):
            resi_dict['class'] = resi.pop()
        if len(resi) > 0:
            if str.isdigit(resi[0]):
                resi_dict['number'] = resi[0]
                del resi[0]
        return resi_dict
    
    @staticmethod
    def get_partnumber(partstring):
        """
        get the part number from a string like PART 1 oder PART 2 -21
        Parameters
        ----------
        partstring: string like 'PART 2 -21'

        >>> FindAtoms.get_partnumber(partstring='PART 2 -21')
        2
        """
        partstring = partstring.upper()
        part = partstring.split()
        try:
            partnum = int(part[1])
        except(ValueError, IndexError):
            print('*** Wrong PART definition found! Check your PART instructions ***')
            partnum = 0
        return partnum
    
    def find_atoms_to_replace(self, frag_atoms, cell, remdist=1.2, only_this=None):
        """
        this method looks around every atom of the fitted fragment and removes
        atoms that are near a certain distance to improve the replace mode

        :param frag_atoms: atoms of the fitting fragment ['C1', '1', 'x', 'y', 'z']
        :param cell: unit cell parameters (list)
        :param remdist: distance below atoms shoud be deleted

        Parameters
        ----------
        only_this: replace only this special atom
        """
        atoms_to_delete = []
        frag_coords = self.get_atomcoordinates(frag_atoms)
        for resinum in self._residues:
            suffix = ''
            if resinum != 0:
                suffix = '_{}'.format(resinum)
            # resinum is the resideue number
            for atom in self._residues[resinum]:
                if atom[0].startswith('Q'):
                    # ignore q peaks:
                    continue
                if only_this and only_this != atom[5]:
                    continue
                # atomlist[4] is the part number
                # replace only in part 0:
                if int(atom[4]) == 0:
                    for name in frag_coords:
                        # name is the atom name
                        if atom[0]+suffix in frag_coords:
                            # do not delete the fitted fragment
                            break
                        at1 = [float(x) for x in frag_coords[name]]
                        at2 = [float(x) for x in atom[1]]
                        resinum1 = self.get_atoms_resinumber(name)
                        resinum2 = self.get_atoms_resinumber(atom[0]+suffix)
                        if at1 == at2 and resinum1 == resinum2:
                            # do not delete atoms on exactly the same position
                            # and same residue
                            break
                        d = atomic_distance(at1, at2, cell)
                        # now get the atom types of the pair atoms and with that
                        # the covalence radius. 
                        if d < remdist:
                            atoms_to_delete.append(atom[0]+suffix)
        return sorted(atoms_to_delete)

    def collect_residues(self):
        """
        find all atoms and sort them into residues

        residues is a dictionary which includes a dictionary for each residue
        which in turn includes a list of its atoms.

        residues = { {'0': ['C1', ['x', 'y', 'z'], linenumber, class, part, element, sfac_number],
                           ['C2', ['x', 'y', 'z'], linenumber, class, part, element, sfac_number]},
                     {'1': ['C1', ['x', 'y', 'z'], linenumber, class, part, element, sfac_number],
                     []} }

        for i in list(residues):
            ats = residues[i]
            for x in ats:
                if 'C12' in x[0]:
                    print x #print C12 from all residues
        """
        resi = False
        part = False
        partnum = '0'
        resiclass = None
        residues = {'0': []}
        resinum = 0
        for num, i in enumerate(self._reslist):
            i = i.upper()
            # First collect the residue
            if re.match(r'^RESI\s+0', i) and resi:
                resi = False
                continue
            if i.startswith(('END', 'HKLF')) and resi:
                resi = False
                continue
            if i.startswith('RESI') and not re.match(r'^RESI\s+0', i):
                resi = True
                resinum = self.get_resi_definition_dict(i.split())['number']
                resiclass = self.get_resi_definition_dict(i.split())['class']
                residues.update({resinum: []})
                continue
            # Now collect the part:
            if re.match(r'^PART\s+0', i) and part:
                part = False
                partnum = 0
                continue
            if i.startswith(('END', 'HKLF')) and part:
                part = False
                continue
            if i.startswith('PART') and not re.match(r'^PART\s+0', i):
                part = True
                partnum = self.get_partnumber(i)
                continue
            #####################
            if resi:
                atom = self.is_atom(i)
                if atom:
                    sfac = atom[1]
                    elem = self.e2s.sfac_2_elem(sfac)
                    residues[resinum].append([atom[0], [float(i) for i in atom[2:5]], num, 
                                              resiclass, partnum, elem, sfac])
            else:
                atom = self.is_atom(i)
                resinum = '0'   # all other atoms are residue 0
                if atom:
                    sfac = atom[1]
                    elem = self.e2s.sfac_2_elem(sfac)
                    residues[resinum].append([atom[0], [float(i) for i in atom[2:5]], num, 
                                              resiclass, partnum, elem, sfac])
        return residues


    def get_atoms_resiclass(self, atom):
        """
        returns the residue class of a given atom. C1 would be 'None'
        C1_1 would be 'CF3' for example
        
        :param atom: an atom name with or without residue number like C1 or C1_1
        :type atom: string
        """
        num = self.get_atoms_resinumber(atom)
        atom = atom.split('_')[0]
        residue = self._residues[num]
        if atom in residue[1]:
            #  atom-name  class
            return residue[1][3]


    def get_atoms_resinumber(self, atom):
        """
        returns the residue number of a given atom. C1 would be '0'
        C1_1 would be '1', ...

        :param atom: an atom name with or without residue number like C1 or C1_1
        :type atom: string
        """
        if '_' in atom:
            suffix = atom.split('_')
            resinum = suffix[-1].strip(string.ascii_letters) # we don't need the part here
            if not resinum:
                resinum = '0'
            if len(resinum) > 4:
                print('*** Invalid residue number in {}'.format(atom))
                sys.exit(-1)
            return str(resinum)
        else:
            return '0'


    def get_atomcoordinates(self, atoms):
        """
        finds an atom regardless if it is in a residue or not.

        rerturns a dictionary {'C1': ['1.123', '0.7456', '3.245']}

        start with digit-> rest auch digit-> resinumber or alias
        start with letter-> rest letter or digit -> residue class

        Negative parts have no letter? part -1 resi 2 is C1_2

        Atom C1 in residue 2 would be addressable with C1_2
        searching for PART is Afaik not neccesary, because inside a residue
        the atom names have to be unique anyway.
        e.g.
        rem DSR put toluene with C1 C2 C3 on C1_2 q2 q3

        :param atoms: list of atoms like ['C1', 'Q2', 'C3_2', ...]
        :type atoms: list
        """
        atom_dict = {}
        for i in atoms:
            num = self.get_atoms_resinumber(i)
            try:
                self._residues[num]
            except KeyError:
                print('*** Atom "{}" not found in res file!! ***'.format(i))
                break
            for x in self._residues[num]:
                if x[0].upper() == i.split('_')[0].upper():
                    single_atom = {i.upper(): x[1]}
                    atom_dict.update(single_atom)
        for i in atoms:
            i = i.upper()
            if i not in list(atom_dict):
                print('\n*** Atom "{}" not found in res file! ***'.format(i))
                sys.exit(0)
        return atom_dict

    def get_atom_line_numbers(self, atoms):
        """
        returns the line numbers in the res_list of a given atom list
        one atom of self._residues[resinum]:
            ['C12', ['0.471727', '0.649578', '0.232054'], 98]
        returns a list like ['98', 'xx', ..]
        :param atoms: list of atom names
        :type atoms: list
        :return lines: list of integers
        """
        # this is to make sure the atomlines are correct: 
        residues = self.collect_residues()
        lines = []
        for at in atoms:
            resinum = self.get_atoms_resinumber(at.upper())
            for x in residues[resinum]:
                atom1 = x[0].upper()
                atom2 = at.split('_')[0].upper()
                if atom1 == atom2:
                    single_atom = x[2] # x[2] is the line number
                    lines.append(single_atom)
        return lines

    def remove_adjacent_hydrogens(self, atoms, sfac_table, atype='H'):
        '''
        if an atom is replaced, its hydrogen atoms are deleted
        this method searches for the first afix behind the atom,
        deletes all lines until AFIX 0, HKLF or the 10th line appears.
        Also Hydrogen atoms outside AFIX are deleted. They are recognized by their
        SFAC number.
        -TODO: create global connectivity table to find atom connections in CF3 CH3 groups
        :param atoms:   list   list of non-Hydrogen atoms
        :param sfac_table: list of atom types ['C', 'H', 'N', ...]
        '''
        delcount = []
        if atype == 'H':
            name = "hydrogen"
        if atype == 'F':
            name = "fluorine"
        lines = self.get_atom_line_numbers(atoms)
        try:
            hydrogen_sfac = sfac_table.index(type)+1
        except(ValueError):
            hydrogen_sfac = False
            return
        for i in lines:
            afix = False
            if i == '':
                continue
            if i == '\n':
                continue
            for n in range(1, 11):
                try:
                    line = self._reslist[i+n+1].upper()
                    atom = line.split()[0].upper()
                    #print('Zeile:', i, 'range:', n, 'line:', line)
                except(IndexError):
                    continue
                if line.startswith('HKLF'):
                    break # stop in this case because the file has ended anyway
                if re.match(atomregex, line) and not afix:
                    # stop if next line is an atom and we are not inside an "AFIX MN"
                    if str(line.split()[1]) == str(hydrogen_sfac):
                        print('Deleted {0} atom {1}'.format(name, atom))
                        delcount.append(atom)
                        self._reslist[i+n] = ''
                        continue
                    else:
                        break
                if line.startswith('AFIX') and afix and line.split()[1] != '0':
                    #print('next afix begins', i)
                    afix = False
                    # stop also if next "AFIX mn" begins
                    break
                if line.startswith('AFIX') and line.split()[1] != '0':
                    #print('AFIX MN starts:', line, 'linenumber:', i+n+1)
                    self._reslist[i+n+1] = ''
                    afix = True # turn on afix flag if first "AFIX mn" is found
                    continue
                if line.startswith('AFIX') and line.split()[1] == '0':
                    #print('AFIX 0:', line)
                    afix = False # turn of afix flag if afix is closed with "AFIX 0"
                    self._reslist[i+n+1] = ''
                    continue
                if afix:
                    try:
                        #print('Atom:', atom)
                        if atom in SHX_CARDS:
                            continue
                        # delete the hydrogen atom
                        self._reslist[i+n+1] = ''
                        print('Deleted {0} atom {1}'.format(name, atom))
                        delcount.append(atom)
                    except(IndexError):
                        continue
        return delcount



def check_source_target(db_source_atoms, res_target_atoms, dbatoms):
    '''
    several checks if the atoms in the dsr command line are consistent
    :param db_source_atoms:   ['C1', 'O1', 'C2', ...]
    :param res_target_atoms:  ['C1', 'Q2', 'C3_2', ...]
    :param dbatoms:           [['C1', 1, '-0.00146', '0.26814', '0.06351'],
                               ['C2', 1, '-1.13341', '-0.23247', '-0.90730'], ...]]
    '''
    temp = [i[0].upper() for i in dbatoms]
    # check if source and target are of same length:
    nsrc = len(db_source_atoms)
    ntrg = len(res_target_atoms)
    if nsrc != ntrg:
        print('*** Number of source and target atoms/peaks is different!! '\
                '({} and {} atoms/peaks) ***'.format(nsrc, ntrg))
        sys.exit(False)
    # do the source atoms exist at all?:
    for i in db_source_atoms:
        i = i.upper()
        if i not in temp:
            print('\n*** Atom {} not found in database entry! ***'.format(i))
            sys.exit(False)
    return True

def rename_dbhead_atoms(new_atoms, old_atoms, dbhead):
    '''
    returns the dbentry header with the old atom names replaced by the new
    number scheme

    :param new_atoms: ['C1A', 'O1A', 'C2A', ...]
    :param old_atoms: ['C1', 'O1', 'C2', ...]
    :param dbhead:    [SAME F5A F6A F6A F4A F7A F8A F8A F9A F9A F7A\n', 'SIMU O1A > F9A\n',
                'RIGU O1A > F9A\n']
    '''
    new = list(reversed(new_atoms))
    headneu = []
    for line in dbhead:
        line = line.split()
        for x, a in enumerate(old_atoms):
            for n,i in enumerate(line):
                if i == a:
                    line[n] = new[x]
        headneu.append(' '.join(line)+'\n')
    return headneu


def set_final_db_sfac_types(db_atom_types, dbatoms, sfac_table):
    '''
    corrects the sfac types of the dbentry according to sfac card of the
    res file
    :param db_atom_types: element names of each atom in the database entry
                          like ['C', 'C', 'N', ... ]
    :type db_atom_types: list
    :param dbatoms: full atoms of the database entry
    :type dbatoms: list
    :param sfac_table: list of scattering factors from SHELXL
    :type sfac_table: list
    '''
    e2s = Elem_2_Sfac(sfac_table)
    atype = list(reversed(db_atom_types))
    for line in dbatoms:                    # go through db entry
        # replace scattering factor (line[1]) with true one
        line[1] = e2s.elem_2_sfac(atype.pop())
    return dbatoms


class SfacTable():
    '''
    Combines the atomtypes of the resfile SFAC and the dbentry to a global
    SFAC table for the new res file.
    '''

    def __init__(self, reslist, dbatom_types):
        '''

        :param reslist:  SHELXL .res file as list
        :param fragment_atom_types:  list ['N', 'C', 'C', 'C']
        :param res_file_name: str file name like 'p21c.res'
        '''
        self._reslist = reslist
        self._db_atom_types = dbatom_types
        self.elements = [x.upper() for x in atoms]


    def set_sfac_table(self):
        '''
        sets the new global sfac table in the res file
        '''
        sfacline = find_multi_lines(self._reslist, r'SFAC\s+[a-zA-Z]+')  # position of the SFAC card
        unitline = find_line(self._reslist, r'UNIT\s+[0-9]+')     # position of the UNIT card
        try:
            if sfacline[-1] > unitline:
                print('*** SFAC in line {} must be defined before UNIT! ***'.format(sfacline[-1]+1))
                sys.exit()
        except():
            pass
        unit = []
        explicit_scat = []
        regular_sfac_line_num = False
        if len(sfacline) == 1:
            # regular SFAC table
            sfac = self._reslist[sfacline[0]].split()[1:]      # SFAC string in the reslist
            sfacline = sfacline[0]
        elif len(sfacline) > 1:
            # first and second type of sfac:
            sfac = []
            for i in sfacline:
                if not ''.join(self._reslist[i].split()).isalpha(): # SFAC with scattering factor
                    # in this case the SFAC command defines also a scattering factor:
                    element = self._reslist[i].split()[1]
                    explicit_scat.append(element) # first SFAC paramter is element
                else:
                    # regular SFAC list
                    sfac.extend(self._reslist[i].split()[1:])
                    regular_sfac_line_num = i
        if regular_sfac_line_num:
            sfacline = regular_sfac_line_num
        if not sfacline:
            print('*** No SFAC card found! Can not proceed ***')
            sys.exit()
        sfac = [x.upper() for x in sfac]
        explicit_scat = [x.upper() for x in explicit_scat]
        for i in self._db_atom_types:  # this is to compare the occurence of element type from resfile and db
            i = i.upper()
            if i not in sfac+explicit_scat:         # all atom types from db not already in sfac
                sfac.append(i)        # get appended to sfac
            if i not in self.elements:
                print('*** Error, atom {} not valid ***'.format(i))
                sys.exit(False)
        for i in range(1, len(sfac+explicit_scat)+1):
            i = str(i).upper()
            unit.append(i)
        # now the sfac and unit tables are written to the resfile
        self._reslist[sfacline] = 'SFAC  {}\n'.format('  '.join(sfac))
        self._reslist[unitline] = 'UNIT  {}\n'.format('  '.join(unit))  # builds the UNIT line
        return sfac+explicit_scat


#############################################################################

class Elem_2_Sfac():
    def __init__(self, sfac_table):
        '''
        SFAC2Element and vice versa conversion
        :param sfac_table: SFAC table of the res file.
        :type sfac_table: list
        '''
        self._sfac_table = sfac_table

    def elem_2_sfac(self, atom_type):
        '''
        returns an sfac-number for the element given in "atom_type"
        :param atom_type: string 'C'
        :type atom_type: string
        '''
        for num, element in enumerate(self._sfac_table, 1):
            if atom_type.upper() == element.upper():
                return num         # return sfac number
                break


    def sfac_2_elem(self, sfacnum):
        '''
        returns an element and needs an sfac-number
        :param sfacnum: string like '2'
        '''
        try:
            sfacnum = int(sfacnum)
        except:
            return False
        for num, element in enumerate(self._sfac_table, 1):
            if sfacnum == num:
                return element.upper()           # return Element name
                break


#############################################################################



class NumberScheme():
    '''
    This class needs a dbentry and the resfile contents to generate a unique
    numbering scheme for the dbentry.
    It essentially adds a letter from the alphabet to the atom name and iterates
    over the names until the names are unique.
    '''

    def __init__(self, reslist, dbatome, resi=False):
        '''
        
        :param reslist: resfile as list
        :param dbatome: list of atoms as list
        :param resi: residue True or false
        '''
        self.__resi = resi
        self._reslist = reslist
        self.__rlist = []
        for i in get_atoms(self._reslist):
            if not re.match('H[0-9]+', i[0]):
                if not re.match('Q[0-9]+', i[0]):
                    self.__rlist.append(i[0])
        self.__dbatome = dbatome
        # Alphabet for the naming scheme:
        self.__alphabet = []
        self.__letters = string.ascii_uppercase
        self.__alphabet = [ i for i in string.ascii_uppercase ]
        # unsorted atomlist:
        self.__atyp = get_atomtypes(self.__dbatome)
        # sorted atomlist:
        self.__atomtype = sorted(self.__atyp[:])


    def _get_number_of_at_in_type(self, atomtype, check):
        '''returns the number of atoms which have the same type'''
        num = 0  # @UnusedVariable
        num = atomtype.count(check)
        return num


    def _generate_numbers(self, atype, num, suffix = ''):
        '''generates numberscheme'''
        atlist = []
        if num == 1:   # in case of only one atom of this type
            atlist.append("%s%s%s"%(atype, num, suffix))
            return atlist
        for i in range(1, num+1):
            atlist.append("%s%s%s"%(atype, str(i), suffix))
        return atlist


    def _contains(self, atomlist, rlist):
        '''checks if one of the atoms in the list exist in resfile'''
        x = [e for e in atomlist if e in '\n'.join(rlist)]
        if x:
            return False
        else:
            return True


    def get_fragment_number_scheme(self, extranames=[]):
        '''
        returns a list of atoms of length len(self.__dbatome) whith a
        naming scheme which fits into the resfile. extraname is a list of atoms 
        which should be also excluded from the new names.
        :param extranames: list of excluding atoms e.g. ['O1A']
        :type extranems: list
        '''
        self.__rlist.extend(extranames)
        if self.__resi:
            print('RESI instruction is enabled. Leaving atom numbers as they are.')
            return [i[0].upper() for i in self.__dbatome]
        num = 0
        atom = []
        newatom = []
        atomlist = []
        for x, i in enumerate(self.__atomtype):
            alph = self.__alphabet[:]
            if atom == i:
                continue
            atom = i
            # number of atoms with given atomtype:
            num = self._get_number_of_at_in_type(self.__atomtype, atom) # unittest: atom instead of atomtype
            # list of atoms with given atomtype
            atomlist = self._generate_numbers(atom, num)
            # important if we want to have all atoms with the same suffix:
            self.__rlist.extend(atomlist)
            # check if rlist contains an atom name of atomlist.
            # If so, add a suffix from aplhabet to it.
            while self._contains(atomlist, self.__rlist) == False:
                suffix = alph[0]
                del alph[0]
                atomlist = self._generate_numbers(atom, num, suffix)
            newatom.extend(atomlist)
        # retain the original order of the atomlist:
        # index of the sorted atoms:
        aindex = sorted(list(range(len(self.__atyp))), key=lambda k: self.__atyp[k])
        orglist = [''] * len(newatom)
        for x, i in enumerate(aindex):
            orglist[i] = newatom[x] # every ith aindex element is replaced by newatom[x] to retain original order
        return orglist


if __name__ == '__main__':
    import doctest
    failed, attempted = doctest.testmod()#verbose=True)
    if failed == 0:
        print('passed all {} tests!'.format(attempted))
    else:
        print('{} of {} tests failed'.format(failed, attempted))
    sys.exit()

    ###########################################################

    from dsrparse import DSR_Parser
    from resfile import ResList
    from dbfile import global_DB
    from resfile import ResListEdit
    res_file = 'p21c.res'
    res_list = ResList(res_file)
    reslist =  res_list.get_res_list()
    find_atoms = FindAtoms(reslist)
    rle = ResListEdit(reslist, find_atoms)
    gdb = global_DB()
    db = gdb.build_db_dict()

    #resiopt = dsr_dict['resi']
    resiopt = False

    fragment = 'oc(cf3)3'
    #fragment = 'PFA1'
    #fragline = gdb.get_fragline_from_fragment(fragment)
    dbatoms = gdb.get_atoms_from_fragment(fragment)
    #dbhead = gdb.get_head_from_fragment(fragment)

    dsrp = DSR_Parser(reslist, rle)
    dsr_dict = dsrp.get_dsr_dict
    num = NumberScheme(reslist, dbatoms, resiopt)
    # das printet auch auf bilschirm:
    numbers = num.get_fragment_number_scheme(extranames=['O1A'])

    print(numbers)

    #  dbtypes = get_atomtypes(dbatoms)


    print('\n')

    #   misc.wrap_headlines(dbhead)

    #  newnames = rename_dbhead_atoms(numbers, dbatoms, dbhead)

    #print(newnames)
    #for i in dbhead:
    #    print(i.strip('\n'))


#    import misc
#    filename = 'testfile.txt'
#    misc.remove_file(filename)
#    try:
#        dfix_file = open(filename, 'w')  # open the ins file
#    except(IOError):
#        print('Unable to write res file!')
#        sys.exit(-1)
#    for i in dbhead:            #modified reslist
#        dfix_file.write("%s" %i)    #write the new file
#    dfix_file.close()


    #sys.exit()

    fa = FindAtoms(reslist)
 #   
  #  atoms = fa.atoms_as_residues
   # for i in atoms:
    #    for y in atoms[i]:
     #       print(y)
    #sys.exit()

    # this might be used to find nearest atoms in same class to make eadp
    # print(fa.get_atoms_resiclass('C1_2'))

    SFAC = ['C', 'H', 'N']

    #sys.exit()

    #  print('Residue dict:', fa.get_resi_definition_dict('RESI 1 TOL'.split()))
    #print( fa.get_atomcoordinates(['Fe1_1', 'C29', 'Q12']) )
    #  print('line number:', fa.get_atom_line_numbers(['C12', 'C333', 'Q12']))
    #
    res_target_atoms = ['c37', 'c5']
    target_lines = fa.get_atom_line_numbers(res_target_atoms)
    for i in target_lines:
        i = int(i)
        rle.remove_line(i, rem=False, remove=False, frontspace=True)
    fa.remove_adjacent_hydrogens(res_target_atoms, SFAC)
    #

    for i in reslist[105:115]:
        print(i.strip('\n'))
    #
    #  print('get_atomtypes', get_atomtypes(dbatoms))
    #  print()
    #  print('get_atomcoordinates:')
    #  coord = fa.get_atomcoordinates(['C1', 'C22', 'C28', 'Q1'])
    #  for i in coord:
    #      print('{}:\t{:<9} {:<9} {:<9}'.format(i, *coord[i]))
    #  print('jhdgfd')
    #
    #  dsrp = DSR_Parser(reslist)
    #  dsr_dict = dsrp.get_dsr_dict

    print(check_source_target(dsr_dict.get('source'), dsr_dict.get('target'), dbatoms))

    #print dbatoms
    num = NumberScheme(reslist, dbatoms, dsr_dict['resi'])
    num.get_fragment_number_scheme()

    dbtypes = get_atomtypes(dbatoms)
    sfac = SfacTable(reslist, dbtypes)
    print('sfac table:')
    print(sfac.set_sfac_table())
    print('######')
    bad_dbatoms = [['lO1', 3, '-0.01453', '1.66590', '1.66590'], ['C1', 1, '-0.00146', '0.26814', '0.06351']]
    print('test')
    #get_atomtypes(bad_dbatoms)



