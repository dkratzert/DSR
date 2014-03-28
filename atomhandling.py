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
from constants import *
from misc import find_line, get_atoms, remove_partsymbol
import textwrap


__metaclass__ = type  # use new-style classes

def get_atomtypes(dbatoms):
    '''
    find all atoms in a list of shelxl format atom lines. 
    returns a list like ['N', 'C', 'C', 'C'].
    '''
    found = []
    # find lines with atoms and see if they are in the atom list
    # print get_atoms(dbatoms)
    for i in dbatoms:
        i = i[0]    # i is the full atom name with number suffix like C1
        atom=''
        for x in i:       # iterate over characters in i
            if re.match("^[A-Za-z]", x):
                atom = atom.upper()+x.upper()  # add characters to atoms until numbers occure
            else:                          # now we have atoms like C, Ca, but also Caaa
                break

        el.atoms = [x.upper() for x in el.atoms]

        if atom[0:2].upper() in el.atoms:    # fixes names like Caaa to be just Ca
            found.append(atom[0:2].upper())  # atoms first, search for all two-letter atoms
            continue
        elif atom[0] in el.atoms:
            found.append(atom[0])  # then for all one-letter atoms
        else:
            print('\n {} is not a valid atom!!\n'.format(atom))
            sys.exit()
    if len(dbatoms) != len(found):    # do we really need this here??
        print("One of the Atoms in the dbentry is not correct! Exiting...")
        sys.exit(-1)
    return found


def check_db_names_for_consistency(dbhead, dbatoms):
    '''Checks if the Atomnames in the restraints of the dbhead are also in
    the list of the atoms of the respective dbentry'''
    
    for i in dbhead:
        if i not in dbatoms:
            print('Atom', i, 'is in a restraint but not in the dbentry!')
        else:
            pass



class FindAtoms():
    '''
    Finds Atoms and creates a data structure
    '''
    
    def __init__(self, reslist):
        self._reslist = reslist
        self._residues = self.collect_residues()

    
    def get_atom(self, atomline):
        '''
        returns all atoms found in the input as list
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
            if re.match(r'^RESI\s+0', i.upper()) and resi:
                resi = False
                continue
            if i.upper().startswith(('END', 'HKLF')) and resi:
                resi = False
                continue
            if i.upper().startswith('RESI') and not re.match(r'^RESI\s+0', i.upper()):
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
            

    def get_atoms_resinumber(self, atom):    
        '''
        returns the residue number of a given atom. C1 would be '0'
        C1_1 would be '1', ...
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
        
        Negative parts habe no letter? part -1 resi 2 is C1_2
        
        Atom C1 in residue 2 would be adressable with C1_2
        searching for PART is Afaik not neccesary, because inside a residue
        the atom names have to be unique anyway.
        e.g.
        rem DSR put toluene with C1 C2 C3 on C1_2 q2 q3
        '''
        atom_dict = {}
        for i in atoms:
            num = self.get_atoms_resinumber(i)
            #print(num, 'rdg')
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
    

    def remove_adjacent_hydrogens(self, atoms):
        '''
        if an atom is replaced, its hydrogen atoms are deleted
        this method searches for the first afix behind the atom, 
        deletes all lines until AFIX 0, HKLF or the 10th line appears.
        '''
        lines = self.get_atom_line_numbers(atoms)
        for i in lines:
            afix = False
            for n in range(0, 10):
                try:
                    line = self._reslist[i+n].upper()
                except(IndexError):
                    continue
                if line.startswith('HKLF'):
                    break # stop in this case because the file has ended anyway
                if line.startswith('AFIX') and afix: # stop also if next afix begins
                    afix = False
                    break
                if line.startswith('AFIX') and line.split()[1] != '0':
                    self._reslist[i+n] = '' 
                    afix = True # turn on afix flag if first afix is found
                    continue
                if line.startswith('AFIX') and line.split()[1] == '0':
                    afix = False # turn of afix flag if afix is closed with "afix 0" 
                    self._reslist[i+n] = ''
                    break
                if afix:
                    try: 
                        atom = line.split()[0].upper()
                        if atom in SHX_CARDS:
                            continue
                        print('Deleted atom {}'.format(line.split()[0]))
                        self._reslist[i+n] = ''
                    except(IndexError):
                        continue
                
    
    
def check_source_target(db_source_atoms, res_target_atoms, dbatoms):
    '''
    several checks if the atoms in the dsr command line are consistent
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
    
  
def wrap_long_lines(head):
    '''
    '''
    pass


def rename_dbhead_atoms(new_atoms, old_atoms, dbhead):
    '''
    returns the dbentry header with the old atom names replaced by the new 
    number scheme
    dbhead = [SAME F5A F6A F6A F4A F7A F8A F8A F9A F9A F7A\n', 'SIMU O1A > F9A\n', 
                'RIGU O1A > F9A\n']
    '''
    new_atoms = list(reversed(new_atoms))
    for x, i in enumerate(old_atoms):
        i = i[0]
        for num, line in enumerate(dbhead):
            line = ' '.join(line.strip().split(' '))
            line = line.replace(i, new_atoms[x])
            dbhead[num] = line+'\n'
    #for num, line in enumerate(dbhead):
    #    line = textwrap.wrap(line, 78, subsequent_indent = '  ')
    #    if len(line) > 1:
    #        line[0] = line[0]+' ='
    #        line[1] = line[1]+'\n'
    #        line = '\n'.join(line)
    #        dbhead[num] = line
    return dbhead



class SfacTable():
    '''
    Combines the atomtypes of the resfile SFAC and the dbentry to a global 
    SFAC table for the new res file.
    '''
    
    def __init__(self, reslist, dbtypes, res_file): 
        self.__reslist = reslist
        self.__db = dbtypes
        self.__resfilename = res_file
    
    
    def set_sfac_table(self):
        '''sets the new global sfac table in the res file'''
        sfacline = find_line(self.__reslist, 'SFAC\s+[a-zA-Z]+')  # position of the SFAC card
        unitline = find_line(self.__reslist, 'UNIT\s+[0-9]+')     # position of the UNIT card
        unit = []
        sfac = self.__reslist[sfacline]      # SFAC string in the reslist
        sfac = sfac.split()
        del sfac[0]
        dbtypes = list(set(self.__db))    # atomtypes in the dbentry

        for i in self.__db:                 # this is to compare the occurence of element type from resfile ant db
            if i not in sfac:             # all atom types from db not already in sfac
                sfac.append(i)        # get appended to sfac
        
        for i in range(len(sfac)):
            i = str(1)        # only unity because we can change this later
            unit.append(i)
        # now the sfac and unit tables are written to the resfile
        self.__reslist[sfacline] = 'SFAC  {}\n'.format('  '.join(sfac))
        self.__reslist[unitline] = 'UNIT  {}\n'.format('  '.join(unit))  # builds the UNIT line
        #print self.__reslist[sfacline]
        #print self.__reslist[unitline]
        return sfac


#############################################################################

class Elem_2_Sfac():
    def __init__(self, sfac):
        self.__sfac = sfac

    def elem_2_sfac(self, atom):
        '''
        returns an sfac-number for the element given in "atom" 
        '''
        for num, element in enumerate(self.__sfac):
            num = num+1
            if atom == element:
                return num         # return sfac number
                break
    
    
    def sfac_2_elem(self, sfacnum):
        '''
        returns an element and needs an sfac-number
        '''
        for num, element in enumerate(self.__sfac):
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
        self.__reslist = reslist
        self.__rlist = []
        for i in get_atoms(self.__reslist):
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
        num = 0
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



#####################################################################################
# This function is not used atm.
def get_next_free_name(res, atomtype):
    '''expects a list of atoms and returns the highest index number
    of an atomtype given as string with atomtype'''
    atomlist = []
    suffix = []
    at = AtomHandling()
    for i in at.get_atoms(res):
        if not re.match('H[0-9]+', i[0]):
            if not re.match('Q[0-9]+', i[0]):
                if i[0].startswith(atomtype):
                    atomlist.append(i)
    # return zero for non-existing atoms. then we can start with the next higher number 1!
    if not atomlist:
        return [0]

    types = at.get_atomtypes(atomlist)
    atoms = list(reversed(atomlist))
    #print types
    #for i in atoms:
    #    print i[0]

    # for every atom, get the full name minus the atom string length:
    for i in types:
        if i == atomtype:
            x = atoms.pop()
            l = len(i)
            suff = x[0][l:].upper()
            suffix.append(suff)


    # this function sorts numbers with strings like C1a, C1b, C2:
    def SortNumericStringList(original):
        newList = [NumericString(x) for x in original]
        newList.sort()
        return [x.rawValue for x in newList]

    sortedlist = SortNumericStringList(suffix)
    return sortedlist[-1]
################################################################################# 


if __name__ == '__main__':
    # for testing:
    from dsrparse import DSR_Parser
    from resfile import ResList
    from dbfile import global_DB
    from resfile import ResList, ResListEdit
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
    
    fragment = 'pfan'
    fragline = gdb.get_fragline_from_fragment(fragment)  
    dbatoms = gdb.get_atoms_from_fragment(fragment)      
    dbhead = gdb.get_head_from_fragment(fragment)
    
    dsrp = DSR_Parser(reslist, rle)
    dsr_dict = dsrp.parse_dsr_line()
    num = NumberScheme(reslist, dbatoms, resiopt)
    num.get_fragment_number_scheme()
    dbtypes = get_atomtypes(dbatoms)
    
    #print(dbhead)
    print('\n')
    
    misc.wrap_headlines(dbhead)
    
    
    
    
    for i in dbhead:
        print(i.strip('\n'))
    
    
    
    
    
    
    
    import misc
    filename = 'testfile.txt'
    misc.remove_file(filename)
    try:
        dfix_file = open(filename, 'w')  # open the ins file
    except(IOError):
        print('Unable to write res file!')
        sys.exit(-1)
    for i in dbhead:            #modified reslist
        dfix_file.write("%s" %i)    #write the new file
    dfix_file.close()
    
    
    sys.exit()
    
    fa = FindAtoms(reslist)
    print('Residue dict:', fa.get_resinum('RESI 1 TOL'.split()))
 #   print fa.get_atomcoordinates(['C12', 'C29', 'Q12'])
    print('line number:', fa.get_atom_line_numbers(['C12', 'C333', 'Q12']))
    
    fa.remove_adjacent_hydrogens(['C2', 'c4', 'c5'])
    
    for num, i in enumerate(reslist):
        if num > 65:
            print(i.strip('\r\n'))
        if num > 85:
            break
    
    print('get_atomtypes', get_atomtypes(dbatoms))
    print()
    print('get_atomcoordinates:')
    coord = fa.get_atomcoordinates(['C1', 'C22', 'C28', 'Q1'])
    for i in coord:
        print('{}:\t{:<9} {:<9} {:<9}'.format(i, *coord[i]))
    print()
    print()
    dsrp = DSR_Parser(reslist)
    dsr_dict = dsrp.parse_dsr_line()
    
    print(check_source_target(dsr_dict.get('source'), dsr_dict.get('target'), dbatoms))
    
    #print dbatoms
    num = NumberScheme(reslist, dbatoms, dsr_dict['resi'])
    num.get_fragment_number_scheme()
    dbtypes = get_atomtypes(dbatoms)
    sfac = SfacTable(reslist, dbtypes)
    print(sfac.set_sfac_table())
