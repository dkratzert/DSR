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
from atoms import Element as el
from misc import find_line, get_atoms
from constants import atomregex, SHX_CARDS
#import textwrap


__metaclass__ = type  # use new-style classes

def get_atomtypes(dbatoms):
    '''
    find all atoms in a list of shelxl format atom lines.
    returns a list like ['N', 'C', 'C', 'C'].

    :param dbatoms:  atoms in a database entry
                     [['O1', 3, '-0.01453', '1.66590', '0.10966'],
                      ['C1', 1, '-0.00146', '0.26814', '0.06351'],
                      [...]
                     ]
    :type dbatoms: list
    '''
    found = []
    # find lines with atoms and see if they are in the atom list
    # print get_atoms(dbatoms)
    elements = [x.upper() for x in el.atoms]
    for i in dbatoms:
        i = i[0].upper()    # i is the full atom name with number suffix like C1
        atom=''
        for x in i:       # iterate over characters in i
            if re.match(r'^[A-Za-z#]', x):
                atom = atom+x      # add characters to atoms until numbers occur
            else:                  # now we have atoms like C, Ca, but also Caaa
                break
        if atom[0:2] in elements:    # fixes names like Caaa to be just Ca
            found.append(atom[0:2])  # atoms first, search for all two-letter atoms
            continue
        elif atom[0] in elements:
            found.append(atom[0])  # then for all one-letter atoms
        else:
            print('\n {} is not a valid atom!!\n'.format(atom))
            sys.exit()
    if len(dbatoms) != len(found):    # do we really need this here??
        print("One of the Atoms in the database entry is not correct! Exiting...")
        sys.exit(-1)
    return found


def check_db_names_for_consistency(dbhead, dbatoms):
    '''Checks if the Atomnames in the restraints of the dbhead are also in
    the list of the atoms of the respective dbentry
    Not used ATM
    '''
    for i in dbhead:
        if i not in dbatoms:
            print('Restraints inconsistent! Atom {} is in a restraint but not in the dbentry!'.format(i))
        else:
            pass



class FindAtoms():
    '''
    Finds Atoms and creates a data structure
    '''

    def __init__(self, reslist):
        '''
        :param reslist: SHELXL .res file as list
        '''
        self._reslist = reslist
        self._residues = self.collect_residues()


    def get_atom(self, atomline):
        '''
        returns all atoms found in the input as list

        :param atomline:  'O1    3    0.120080   0.336659   0.494426  11.00000   0.01445 ...'
        '''
        atom = ''
        if re.search(atomregex, str(atomline)):        # search atoms
            atom = atomline.split()[:5]              # convert to list and use only first 5 columns
            if atom[0].upper() not in SHX_CARDS:      # exclude all non-atom cards
                return atom
            else:
                return


    def get_resinum(self, resi):
        '''
        returns the residue number and class of a string like 'RESI TOL 1'
        or 'RESI 1 TOL'
        {'class': 'TOL', 'number': '1'}

        :param resi: 'RESI number class'
        '''
        resi_dict = {
            'class' : None,
            'number': None}
        resi.remove('RESI')
        resi.sort()
        if str.isalpha(resi[-1][0]):
            resi_dict['class'] = resi.pop()
        if len(resi) > 0:
            if str.isdigit(resi[0]):
                resi_dict['number'] = resi[0]
                del resi[0]
        return resi_dict


    def collect_residues(self):
        '''
        find all atoms and sort them into residues

        residues is a dictionary which includes a dictionary for each residue
        which in turn includes a list of its atoms.

        residues = { {'0': ['C1', 'x y z', 'linenumber'], ['C2', 'x y z', 'linenumber']},
                     {'1': ['C1', 'x y z', 'linenumber'], []} }

        for i in residues.keys():
            ats = residues[i]
            for x in ats:
                if 'C12' in x[0]:
                    print x #print C12 from all residues
        '''
        resi = False
        resiclass = None
        residues = {'0': []}
        for num, i in enumerate(self._reslist):
            i = i.upper()
            if re.match(r'^RESI\s+0', i) and resi:
                resi = False
                continue
            if i.startswith(('END', 'HKLF')) and resi:
                resi = False
                continue
            if i.startswith('RESI') and not re.match(r'^RESI\s+0', i):
                resi = True
                resinum = self.get_resinum(i.split())['number']
                resiclass = self.get_resinum(i.split())['class']
                residues.update({resinum: []})
                continue
            if resi:
                atom = self.get_atom(i)
                if atom:
                    residues[resinum].append([atom[0], atom[2:5], num, resiclass])
            else:
                atom = self.get_atom(i)
                resinum = '0'   # all other atoms are residue 0
                if atom:
                    residues[resinum].append([atom[0], atom[2:5], num, resiclass])
        return residues


    def get_atoms_resiclass(self, atom):
        '''
        returns the residue class of a given atom. C1 would be 'None'
        C1_1 would be 'CF3' for example

        :param atom: an atom name with or without residue number like C1 or C1_1
        '''
        #all_resi_numbers = self._residues.keys()
        num = self.get_atoms_resinumber(atom)
        atom = atom.split('_')[0]
        #print(all_resi_numbers)
        residue = self._residues[num]
        if atom in residue[1]:
            #  atom-name  class
            return residue[1][3]


    def get_atoms_resinumber(self, atom):
        '''
        returns the residue number of a given atom. C1 would be '0'
        C1_1 would be '1', ...

        :param atom: an atom name with or without residue number like C1 or C1_1
        '''
        if '_' in atom:
            suffix = atom.split('_')
            resinum = suffix[-1].strip(string.ascii_letters) # we don't need the part here
            if not resinum:
                resinum = '0'
            if len(resinum) > 4:
                print('Invalid residue number in', atom)
                sys.exit(-1)
            return str(resinum)
        else:
            return '0'


    def get_atomcoordinates(self, atoms):
        '''
        finds an atom regardless if it is in a residue or not.

        start with digit-> rest auch digit-> resinumber or alias
        start with letter-> rest letter or digit -> residue class

        Negative parts have no letter? part -1 resi 2 is C1_2

        Atom C1 in residue 2 would be addressable with C1_2
        searching for PART is Afaik not neccesary, because inside a residue
        the atom names have to be unique anyway.
        e.g.
        rem DSR put toluene with C1 C2 C3 on C1_2 q2 q3

        :param atoms: list of atoms like ['C1', 'Q2', 'C3_2', ...]
        '''
        atom_dict = {}
        for i in atoms:
            num = self.get_atoms_resinumber(i)
            try:
                self._residues[num]
            except(KeyError):
                print('Target atom "{}" not found in res file!!'.format(i))
                break
                #return
            for x in self._residues[num]:
                if x[0].upper() == i.split('_')[0].upper():
                    single_atom = {i.upper(): x[1]}
                    atom_dict.update(single_atom)
        for i in atoms:
            i = i.upper()
            if i not in list(atom_dict.keys()):
                print('\nTarget atom "{}" not found in res file!'.format(i))
                sys.exit()
        return atom_dict


    def get_atom_line_numbers(self, atoms):
        '''
        returns the line numbers in the res_list of a given atom list
        one atom of self._residues[resinum]:
            ['C12', ['0.471727', '0.649578', '0.232054'], 98]
        '''
        lines = []
        for at in atoms:
            resinum = self.get_atoms_resinumber(at.upper())
            for x in self._residues[resinum]:
                if x[0].upper() == at.split('_')[0].upper():
                    single_atom = x[2] # x[2} is the line number
                    lines.append(single_atom)
        return lines


    def remove_adjacent_hydrogens(self, atoms, sfac_table):
        '''
        if an atom is replaced, its hydrogen atoms are deleted
        this method searches for the first afix behind the atom,
        deletes all lines until AFIX 0, HKLF or the 10th line appears.
        Also Hydrogen atoms outside AFIX are deleted. They are recognized by their
        SFAC number.
        :param atoms:   list   list of non-Hydrogen atoms
        :param sfac_table: list of atom types ['C', 'H', 'N', ...]
        '''
        lines = self.get_atom_line_numbers(atoms)
        try:
            hydrogen_sfac = sfac_table.index('H')+1
        except(ValueError):
            hydrogen_sfac = False
            return

        for i in lines:
            afix = False
            if i == '':
                continue
            if i == '\n':
                continue
            for n in range(0, 10):
                try:
                    line = self._reslist[i+n].upper()
                    atom = line.split()[0].upper()
                    #print('Zeile:', i, 'range:', n, 'line:', line)
                except(IndexError):
                    continue
                if line.startswith('HKLF'):
                    break # stop in this case because the file has ended anyway
                if re.match(atomregex, line) and not afix:
                    # stop if next line is an atom and we are not inside an "AFIX MN"
                    if str(line.split()[1]) == str(hydrogen_sfac):
                        print('Deleted hydrogen atom {}'.format(atom))
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
                    #print('AFIX MN starts:', line)
                    self._reslist[i+n] = ''
                    afix = True # turn on afix flag if first "AFIX mn" is found
                    continue
                if line.startswith('AFIX') and line.split()[1] == '0':
                    #print('AFIX 0:', line)
                    afix = False # turn of afix flag if afix is closed with "AFIX 0"
                    self._reslist[i+n] = ''
                    continue
                if afix:
                    try:
                        #print('Atom:', atom)
                        if atom in SHX_CARDS:
                            continue
                        # delete the hydrogen atom
                        self._reslist[i+n] = ''
                        print('Deleted Hydrogen atom {}'.format(atom))
                    except(IndexError):
                        continue



def check_source_target(db_source_atoms, res_target_atoms, dbatoms):
    '''
    several checks if the atoms in the dsr command line are consistent
    :param db_source_atoms:   ['C1', 'O1', 'C2', ...]
    :param res_target_atoms:  ['C1', 'Q2', 'C3_2', ...]
    :param dbatoms:           ['C1', 'C2', 'C3', ...]
    '''
    dbat = []
    for i in dbatoms:
        dbat.append(i[0].upper())

    # check if source and target are of same length:
    nsrc = len(db_source_atoms)
    ntrg = len(res_target_atoms)
    if nsrc != ntrg:
        print('Number of source and target atoms is different!! '\
                '({} and {})'.format(nsrc, ntrg))
        sys.exit()

    # do the source atoms exist at all?:
    for i in db_source_atoms:
        if i.upper() not in dbat:
            print('\nAtom {} not found in database entry! Exiting...\n'.format(i))
            sys.exit(-1)

    return 'check_source_target() succeded!\n'





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



class SfacTable():
    '''
    Combines the atomtypes of the resfile SFAC and the dbentry to a global
    SFAC table for the new res file.
    '''

    def __init__(self, reslist, fragment_atom_types, res_file_name):
        '''

        :param reslist:  SHELXL .res file as list
        :param fragment_atom_types:  list ['N', 'C', 'C', 'C']
        :param res_file_name: str file name like 'p21c.res'
        '''
        self._reslist = reslist
        self._db_atom_types = fragment_atom_types
        self.__res_file_name = res_file_name


    def set_sfac_table(self):
        '''
        sets the new global sfac table in the res file
        '''
        sfacline = find_line(self._reslist, r'SFAC\s+[a-zA-Z]+')  # position of the SFAC card
        unitline = find_line(self._reslist, r'UNIT\s+[0-9]+')     # position of the UNIT card
        unit = []
        sfac = self._reslist[sfacline]      # SFAC string in the reslist
        sfac = sfac.split()
        del sfac[0]
        #dbtypes = list(set(self._db_atom_types))    # atomtypes in the dbentry

        for i in self._db_atom_types:                 # this is to compare the occurence of element type from resfile ant db
            if i not in sfac:             # all atom types from db not already in sfac
                sfac.append(i)        # get appended to sfac

        for i in range(len(sfac)):
            i = str(1)        # only unity because we can change this later
            unit.append(i)
        # now the sfac and unit tables are written to the resfile
        self._reslist[sfacline] = 'SFAC  {}\n'.format('  '.join(sfac))
        self._reslist[unitline] = 'UNIT  {}\n'.format('  '.join(unit))  # builds the UNIT line
        #print self._reslist[sfacline]
        #print self._reslist[unitline]
        return sfac


#############################################################################

class Elem_2_Sfac():
    def __init__(self, sfac):
        self._sfac = sfac

    def elem_2_sfac(self, atom):
        '''
        returns an sfac-number for the element given in "atom"
        :param atom: string 'C1'
        '''
        for num, element in enumerate(self._sfac):
            num = num+1
            if atom == element:
                return num         # return sfac number
                break


    def sfac_2_elem(self, sfacnum):
        '''
        returns an element and needs an sfac-number
        :param sfacnum: string like '2'
        '''
        for num, element in enumerate(self._sfac):
            num = num+1
            if sfacnum == num:
                return element           # return Element name
                break


#############################################################################



class NumberScheme():
    '''
    This class needs a dbentry and the resfile contents to generate a unique
    numbering scheme for the dbentry.
    It essentially adds a letter from the alphabet to the atom name and iterates
    over the names until the names are unique.
    '''

    def __init__(self, reslist, dbatome, resi):
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


    def get_number_of_at_in_type(self, atomtype, check):
        '''returns the number of atoms which have the same type'''
        num = 0  # @UnusedVariable
        num = atomtype.count(check)
        return num


    def generate_numbers(self, atype, num, suffix = ''):
        '''generates numberscheme'''
        atlist = []
        if num == 1:   # in case of only one atom of this type
            atlist.append("%s%s%s"%(atype, num, suffix))
            return atlist
        for i in range(1, num+1):
            atlist.append("%s%s%s"%(atype, str(i), suffix))
        return atlist


    def contains(self, atomlist, rlist):
        '''checks if one of the atoms in the list exist in resfile'''
        x = [e for e in atomlist if e in '\n'.join(rlist)]
        if x:
            return False
        else:
            return True


    def get_fragment_number_scheme(self):
        ''' returns a list of atoms of length len(self.__dbatome) whith a
            naming scheme which fits into the resfile.
        '''
        if self.__resi:
            print('RESI instruction is enabled. Leaving atom numbers as they are.')
            atoms = []
            for i in self.__dbatome:
                atoms.append(i[0].upper())
            return atoms
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
            num = self.get_number_of_at_in_type(self.__atomtype, atom) # unittest: atom instead of atomtype
            # list of atoms with given atomtype
            atomlist = self.generate_numbers(atom, num)
            # important if we want to have all atoms with the same suffix:
            self.__rlist.extend(atomlist)
            # check if rlist contains an atom name of atomlist.
            # If so, add a suffix from aplhabet to it.
            while self.contains(atomlist, self.__rlist) == False:
                suffix = alph[0]
                del alph[0]
                atomlist = self.generate_numbers(atom, num, suffix)
            newatom.extend(atomlist)

        # retain the original order of the atomlist:
        # index of the sorted atoms:
        aindex = sorted(list(range(len(self.__atyp))), key=lambda k: self.__atyp[k])
        orglist = [''] * len(newatom)
        for x, i in enumerate(aindex):
            orglist[i] = newatom[x] # every ith aindex element is replaced by newatom[x] to retain original order

        print('Fragment atom names:', ', '.join(orglist))
        return orglist



if __name__ == '__main__':
    # for testing:
    from dsrparse import DSR_Parser
    from resfile import ResList
    from dbfile import global_DB
    from resfile import ResListEdit
    import misc
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
    fragline = gdb.get_fragline_from_fragment(fragment)
    dbatoms = gdb.get_atoms_from_fragment(fragment)
    dbhead = gdb.get_head_from_fragment(fragment)

    dsrp = DSR_Parser(reslist, rle)
    dsr_dict = dsrp.parse_dsr_line()
    num = NumberScheme(reslist, dbatoms, resiopt)
    # das printet auch auf bilschirm:
    numbers = num.get_fragment_number_scheme()
    #print('#######################')
    print(numbers)
    #print('#######################')
    dbtypes = get_atomtypes(dbatoms)


    print('\n')

    misc.wrap_headlines(dbhead)

    newnames = rename_dbhead_atoms(numbers, dbatoms, dbhead)

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

    # this might be used to find nearest atoms in same class to make eadp
    print(fa.get_atoms_resiclass('C1_2'))

    SFAC = ['C', 'H', 'N']

    #sys.exit()

  #  print('Residue dict:', fa.get_resinum('RESI 1 TOL'.split()))
  #  print fa.get_atomcoordinates(['C12', 'C29', 'Q12'])
  #  print('line number:', fa.get_atom_line_numbers(['C12', 'C333', 'Q12']))
  #
  #  fa.remove_adjacent_hydrogens(['C2', 'c4', 'c5'], SFAC)
  #
  # # for num, i in enumerate(reslist):
  # #     if num > 65:
  # #         print(i.strip('\r\n'))
  # #     if num > 85:
  # #         break
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
  #  dsr_dict = dsrp.parse_dsr_line()

    print(check_source_target(dsr_dict.get('source'), dsr_dict.get('target'), dbatoms))

    #print dbatoms
    num = NumberScheme(reslist, dbatoms, dsr_dict['resi'])
    num.get_fragment_number_scheme()
    dbtypes = get_atomtypes(dbatoms)
    print('#############', dbtypes, '########dbtypes###################')
   # sfac = SfacTable(reslist, dbtypes)
   # print(sfac.set_sfac_table())
