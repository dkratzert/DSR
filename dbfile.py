# -*- encoding: utf-8 -*-
# möp
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

from collections import Counter
import os, sys
import re
import tarfile

from constants import atomregex, SHX_CARDS, RESTRAINT_CARDS, sep_line
from atoms import Element
from atomhandling import get_atomtypes
from misc import atomic_distance, nalimov_test, std_dev, median, pairwise,\
    unwrap_head_lines
from copy import deepcopy
from os.path import expanduser
from misc import touch



def invert_dbatoms_coordinates(atoms):
    '''
    Inverts SHELXL atom coordinates like
    [[C1  1  0.44  0.21  -1.23 ][ ...]]

    :param atoms: list of atom list
    '''
    for line in atoms:
        try:
            inv_coord = [ str(-float(i)) for i in line[-3:] ]
        except:
            print('Unable to invert fragment coordinates.')
            return False
        line[-3:] = inv_coord
    return atoms


def search_fragment_name(search_string):
    '''
    searches the Name: comments in the database for a given name
    '''
    from misc import dice_coefficient
    gdb = global_DB()
    db = gdb.build_db_dict()
    frags = db.keys()
    names_list = []
    for i in frags:
        fragname = gdb.get_comment_from_fragment(i)
        line_number = gdb.get_line_number_from_fragment(i)
        dbname = gdb.get_db_name_from_fragment(i)
        names_list.append([i, fragname, line_number, dbname])
    search_results = {}
    for i in names_list:
        db_entry = i[1]
        #Levenshtein gibt bei kurzen Suchstrings zu schlechte Ergebnisse:
        #coefficient = levenshtein(self.search_string, db_entry)
        coefficient = dice_coefficient(search_string, db_entry)
        search_results[coefficient] = i
    # select the best 5 results:
    selected_results = [search_results[i] for i in sorted(search_results)[0:5]]
    #for i in selected_results:
    #    i.append(make_sortkey(i[1]))
    #selected_results.sort(key=lambda x: x[4].lower())
    return selected_results


def print_search_results(results):
    '''
    prints the results of a database search to screen and exit.
    results are
    '''
    #print(' Found following database entries:')
    print(' Fragment          | Full name, Comments                      | Line number')
    print(sep_line)
    for line in results:
        print(' {:15s}   | {:40s} | {}'.format(line[0], line[1], line[2]))
    sys.exit()

def make_sortkey(full_name):
    """
    Algorythm inspired by W. Sage J. Chem. Inf: Comput. Sci. 1983, 23, 186-197
    """
    full_name = ''.join(e for e in full_name if e not in ('{}()[]'))
    full_name = full_name.split(',')[0]
    if full_name.startswith('tert-'):
        full_name = full_name[4:]
    if full_name.startswith('sec-'):
        full_name = full_name[3:]
    if full_name.startswith('t-'):
        full_name = full_name[1:]
    if full_name.startswith('p-'):
        full_name = full_name[2:]
    if full_name.startswith('o-'):
        full_name = full_name[1:]
    if full_name.startswith('n-'):
        full_name = full_name[1:]
    if full_name.startswith('m-'):
        full_name = full_name[1:]
    if full_name.startswith('iso-'):
        full_name = full_name[3:]
    if full_name.startswith('i-'):
        full_name = full_name[1:]
    full_name = ''.join(e for e in full_name if e not in ('-_.\'1234567890, '))
    return full_name

__metaclass__ = type  # use new-style classes

# hardwired names of the database files:
# dsr_db.txt is the db from the distribution. This file should not be edited.
# dsr_user_db.txt if this file exists, all its content is also read in.

class ReadDB():
    '''
    reads in the system and user db files (self._db_file_names) and makes
    a dictionary of them.
    '''

    def __init__(self, main_dbdir='./'):
        self.maindb = "dsr_db.txt" 
        self.userdb = "dsr_user_db.txt"
        self.homedir = expanduser("~")
        self._db_dir = main_dbdir
        self._databases = self.getDB_files_dict(self.maindb, self.userdb)

    @property
    def get_databases(self):
        return self._databases

    def getDBpath(self, db_dir, db_file_name):
        '''
        returns the full db path as tuple
        '''
        fullpath = os.path.join(db_dir, db_file_name)  # full path with filename
        return fullpath


    def getDB_files_dict(self, maindb='', userdb=''):
        '''
        returns the database as dictionary. Each file has its own key.
        {'dsr-db': ('line1\n', 'line2\n', '...'), 'dsr-user-db': ('line1\n', 'line2\n', '...')}
        '''
        db_dict = {}
        for name in [maindb, userdb]:
            dblist = [] 
            if name == userdb:
                filepath = self.getDBpath(self.homedir, name)
                if not os.path.isfile(filepath):
                    touch(filepath)
            else:
                filepath = self.getDBpath(self._db_dir, name)
            base_filename = os.path.splitext(name)[0]
            try:
                with open(filepath, 'r') as f:
                    for line in f:
                        if line.startswith('#'):
                            line = ''
                        dblist.append(line)
            except(IOError) as e:
                print(e)
                sys.exit(-1)
            db_dict[base_filename] = dblist
            del dblist
        return db_dict


    def find_db_tags(self):
        '''
        This method returns all fragment name tags in the database
        '''
        regex = r'^<[^/].*>'  # regular expression for db tag.
        dbnames = []
        for db in self._databases:
            for num, line in enumerate(self._databases[db]):
                if re.match(regex, line):
                    frag = [str(line.strip('<> \n\r')).upper(), num + 1, db]  # stripping spaces is important here!
                    if frag[0] in ['CF3', 'CF6', 'CF9']:
                        print('\nWarning, {} is a reserved fragment name. Ignoring database entry {}.'.format(frag[0], frag[0]))
                        continue
                    dbnames.append(frag)
        nameset = []
        for i in dbnames:
            nameset.append(i[0])
        if len(set(nameset)) != len(nameset):
            c1 = Counter(nameset)
            c2 = Counter(set(nameset))
            diff = c1 - c2
            duplicates = list(diff.elements())
            for i in duplicates:
                # print list(set([ i for y in dbnames if y for y in i]))[0]
                # print dbnames
                print('\nDuplicate database entry "{}" found! Please remove/rename '\
                    'second entry\nand/or check all end tags in the database dsr_usr_db.txt or dsr_db.txt.\n'.format(duplicates.pop()))
            sys.exit(False)
        # # sort lower-case:
        dbnames.sort(key=lambda x: x[0].lower())
        return dbnames



class global_DB():
    '''
    creates a final dictionary where all dbs are included
    '''
    def __init__(self, invert=False, fragment=None):
        '''
        self._db_tags: ['12-DICHLOROBENZ', 590, 'dsr_db']
        
        self._db_plain_dict: dictionary with plain text of the individual databases
                  {'dsr_db': ['# Fragment database of DSR.py\n', '#\n',
                              '# Some Fragments are from geometry optimizations with Gaussian 03:\n',
                              '#  Gaussian 03, Revision B.04,\n', '#  M. J. Frisch, et. al. Gaussian, ...}
        
        self._dbentry_dict: dictionary with the individial fragments
                  {'benzene':
                    {'comment': ['Source: GRADE import', 'Name: Benzene, C6H6'],
                     'head': ['DFIX 1.379 0.015 C1 C2 ', 'DFIX 1.379 0.015 C1 C6 ',
                             'DFIX 1.379 0.015 C2 C3 ', 'DFIX 1.379 0.015 C3 C4 ',
                             'DFIX 1.379 0.015 C4 C5 ', 'DFIX 1.379 0.015 C5 C6 ',
                             'DANG 2.386 0.022 C2 C6 ', 'DANG 2.386 0.022 C1 C3 ',
                             'DANG 2.386 0.022 C2 C4 ', 'DANG 2.386 0.022 C3 C5 ',
                             'DANG 2.386 0.022 C4 C6 ', 'DANG 2.386 0.022 C1 C5 ',
                             'FLAT C1 C2 C3 C4 ', 'FLAT C2 C3 C4 C5 ', 'FLAT C3 C4 C5 C6 ',
                             'FLAT C4 C5 C6 C1 ', 'FLAT C5 C6 C1 C2 ', 'FLAT C6 C1 C2 C3'],
                     'resi': 'BEN1',
                     'name': 'BENZENE',
                     'line': 762,
                     'db': 'dsr_db',
                     'fragline': ['FRAG', '17', '1', '1', '1', '90', '90', '90'],
                     'atoms': [['C1', '1', '-0.012', '-0.0', '-0.007'],
                               ['C2', '1', '-0.012', '-0.0', '-1.388'],
                               ['C3', '1', '-1.208', '-0.0', '-2.078'],
                               ['C4', '1', '-2.404', '-0.0', '-1.388'],
                               ['C5', '1', '-2.404', '-0.0', '-0.007'],
                               ['C6', '1', '-1.208', '-0.0', '0.683']]},
                    'bipyridine':
                     {'comment': ...
        :param invert: inverts the coordinates of a fragment
        :type invert:  boolean
        :param main_dbdir:  directory where the database is located. Default is environment variable DSR_DB_DIR
        :type main_dbdir:   string
        :param dbnames: file names of the databases
        :type dbnames:  string
        '''
        self.invert = invert
        try:
            main_dbdir = os.environ["DSR_DIR"]
        except(KeyError):
            main_dbdir='./'
        self._getdb = ReadDB(main_dbdir)
        self._db_tags = self._getdb.find_db_tags()
        if fragment:
            for num, i in enumerate(self._db_tags):
                if i[0].lower() == fragment.lower():
                    self._db_tags = [self._db_tags[num]] # speedup in case the fragment is known
                    break
        self._db_plain_dict = self._getdb.get_databases
        self._dbentry_dict = self.build_db_dict()

    def list_fragments(self):
        '''
        list all fragments in the db as list of lists
        [['tbu-c', 1723, 'dsr_db', 'Tert-butyl-C'], ...]
        '''
        fraglist = []
        fragments = self._dbentry_dict.keys()
        for frag in fragments:
            comment = self.get_comment_from_fragment(frag)
            line = [frag, 
                    self.get_line_number_from_fragment(frag), 
                    self.get_db_name_from_fragment(frag), 
                    comment, 
                    make_sortkey(comment)]
            fraglist.append(line)
        fraglist.sort(key=lambda x: x[4].lower())
        return fraglist 
    

    def build_db_dict(self):
        '''
        returns the global db dictionary
        '''
        db_dict = {}
        for i in self._db_tags:
            fragment = i[0].lower()
            line = i[1]   # line number where dbentry is in dbfile
            db = str(i[2])
            headboth = self.get_head_lines(fragment, db, line)
            header = headboth[0]
            fragline = headboth[1]
            residue = self.get_residue_from_head(header, fragment)
            atoms = self.get_fragment_atoms(fragment, db, line)
            comment = headboth[2]
            comment = [' '.join(x) for x in comment]
            db_dict[fragment] = {
                'head'    : header, # header with just the restraints
                'resi'    : residue, # the residue class
                'fragline': fragline.split(), # FRAG ...
                'atoms'   : atoms, # the atoms as lists of list
                'line'    : line, # line number
                'db'      : db,   # user or dsr db
                'comment' : comment, # the comment line
                'name'    : i[0]
                }
        if not db_dict:
            print('No database dsr_db.txt found!\n')
            sys.exit()
        return db_dict
    
    @property
    def db_dict(self):
        return self._dbentry_dict

    def get_residue_from_head(self, head, fragment=''):
        '''
        extracts the residue from the head
        and returns just the first string after RESI
        :param head: ['line1\n', 'line2\n', '...']
        :type head: list of strings
        :return residue class
        '''
        for index, line in enumerate(head):
            line = line.upper()
            if line.startswith('RESI'):
                resiline = line.split()
                del head[index]  # remove resi from head
                for n, i in enumerate(resiline):
                    if not i[0].isalpha():
                        del resiline[n]
                try:
                    return resiline[1].upper()
                except(IndexError):
                    print('Invalid residue definition in database entry {}.'.format(fragment))
                    sys.exit(False)
        return False


    def get_fragment_atoms(self, fragment, db_name, line):
        '''
        returns the atoms of a fragment as list
        [['C1', '1', '7.600', '-1.044', '4.188'], [...]]
        :param fragment: db fragment name
        :type fragment:  string
        :param db_name:  database from the database dictionary
        :type db_name:   string
        :param line:     line number of db entry
        :type line:      integer
        '''
        atoms = []
        end = False
        regex = re.escape(r'</{}>'.format(fragment.lower()))
        for i in self._db_plain_dict[db_name][int(line):]:
            i = i.strip('\n\r')
            if re.match(regex, i.lower()):  # find the endtag of db entry
                end = True
                break
            if re.search(atomregex, str(i)):  # search atoms
                l = i.split()[:5]  # convert to list and use only first 5 columns
                if l[0].upper() not in SHX_CARDS:  # exclude all non-atom cards
                    atoms.append(l)
        if not atoms:
            print('Database entry of "{}" in line {} of "{}.txt" is corrupt. No atoms found!'.format(fragment, line, db_name))
            print('Have you really followed the syntax?')
            sys.exit(False)
        if not end:
            print('Could not find end of dbentry for fragment '\
                '"{}" in line {} of "{}.txt". Check your database.'.format(fragment, line, db_name))
            sys.exit(False)
        if self.invert:
            atoms = invert_dbatoms_coordinates(atoms)
        return atoms

    def get_head_for_gui(self, fragment):
        '''
        returns header information of the specific fragment:
        tag, Name/comment, source, cell, residue, dbtype, restr, atoms
        ----------------------------------------
        <tag>
         toluene 
        </tag>
        <comment>
         Toluene, C7H8 
        </comment>
        <source>
         CCDC CESLUJ 
        </source>
        <cell>
         1;;1;;1;;90;;90;;90 
        </cell>
        <residue>
         TOL 
        </residue>
        <dbtype>
         dsr_db 
        </dbtype>
        <restr>
         SADI C2 C3 C3 C4 C4 C5 C5 C6 C6 C7 C7 C2;;SADI 0.04 C2 C6 C2 C4 C7 C5 C3 C7 C4 C6 C3 C5;;DFIX 1.51 C1 C2;;SADI 0.04 C1 C7 C1 C3;;FLAT C1 > C7;;SIMU C1 > C7;;RIGU C1 > C7 
         </restr>
        <atoms>
         C1;;6;;1.78099;;7.14907;;12.00423
        C2;;6;;2.20089;;8.30676;;11.13758
        C3;;6;;1.26895;;9.02168;;10.39032
        C4;;6;;1.64225;;10.07768;;9.58845
        C5;;6;;2.98081;;10.44432;;9.51725
        C6;;6;;3.92045;;9.74974;;10.25408
        C7;;6;;3.53891;;8.69091;;11.05301 
        </atoms>
        '''
        fragment = fragment.lower()
        print("<tag>\n", fragment, "\n</tag>")
        print("<comment>\n", self.get_comment_from_fragment(fragment), "\n</comment>")
        print("<source>\n", self.get_src_from_fragment(fragment), "\n</source>")
        # it returns cartesian coordinates, so we need this cell here:
        print("<cell>\n", '1;;1;;1;;90;;90;;90', "\n</cell>") 
        print("<residue>\n", self.get_resi_from_fragment(fragment), "\n</residue>")
        print("<dbtype>\n", self.get_db_name_from_fragment(fragment), "\n</dbtype>")
        print('<restr>\n', ';;'.join(self.db_dict[fragment]['head']), '\n', '</restr>')
        self.check_consistency(fragment) # too many critical errors with GUI
        self.check_db_header_consistency(fragment)
        if not self.check_sadi_consistence(fragment):
            sys.exit()
        self.check_db_atom_consistency(fragment)    

    def check_consistency(self, fragment):
        '''
        check if the fragline makes sense and if the fragment_dict is complete
        '''
        try:
            dbentry = self.db_dict[fragment.lower()]
        except KeyError:
            print("Fragment {} not found in the database.".format(fragment))
            sys.exit()
        for i in dbentry:
            if i == 'comment':
                continue
            if not dbentry[i]:
                if i == 'head':
                    print('- Restraints in the header of database entry "{}" missing!\n Check your Database!'\
                        ''.format(fragment))
                    return False
                else:
                    print('- Values for "{}" in database entry "{}" missing! Check your Database!'\
                        ''.format(str.upper(i), fragment))
                    return False
        if len(dbentry['fragline']) != 8:
            print('- The line starting with "FRAG" in the database entry of {} is '\
                'not correct.\n  Are the cell parameters really correct? '\
                '"FRAG 17 a b c alpha beta gamma"\n'.format(fragment))
            sys.exit(False)
        return True

    def get_sum_formula(self, fragment):
        '''
        returns the sum formula of a fragment
        as string like 'SI4 C9 O1'
        :param fragment: fragment name
        :type fragment: string
        '''
        try:
            types = get_atomtypes(self.db_dict[fragment]['atoms'])
        except:
            return None
        formula = ''
        for el in set(types):
            num = types.count(el)
            formula+='{}{} '.format(el, num)
        return formula

    def check_db_atom_consistency(self, fragment):
        '''This method is for atoms only (without db header)!
          check the db for duplicates:
        '''
        dbatoms = [i[0] for i in self.get_atoms_from_fragment(fragment)]
        # check for duplicates:
        while dbatoms:
            at = dbatoms.pop()
            if at in dbatoms:
                print('duplicate atom {0} in database entry "{1}" '\
                    'found!'.format(at, fragment))
                print("Check your database...")
                sys.exit(-1)

    def check_db_header_consistency(self, fragment):
        '''
        - Checks if the Atomnames in the restraints of the dbhead are also in
          the list of the atoms of the respective dbentry.
        - Checks wether restraints cards are vaid.
        '''
        fragment = fragment.lower()
        restraints = self.db_dict[fragment]['head']
        atoms = self.get_atoms_from_fragment(fragment)
        db = self.get_db_name_from_fragment(fragment)
        status = True
        atoms = [i[0].upper() for i in atoms]
        restraint_atoms_list = set([])
        for line in restraints:
            if not line:
                continue
            line = line.upper()
            line2 = line.split()
            # only the first 4 characters, because SADI_TOL would be bad:
            if line2[0] not in SHX_CARDS:  
                status = False
                print('Bad line in header of database entry "{}" found! ({}.txt)'.format(fragment, db))
                print(line)
                #sys.exit(status)
            if line[:4] in RESTRAINT_CARDS:
                line = line[5:].split()
                for i in line:
                    if i in ('>', '<'):
                        continue
                    try:
                        float(i)
                    except(ValueError):
                        restraint_atoms_list.add(i)
        for atom in restraint_atoms_list:
            atom = atom.upper()
            if not atom in atoms:
                status = False
                print('\nUnknown atom "{}" in restraints of "{}".'.format(atom, fragment))
        if not status:
            print('Check database entry.\n')
            sys.exit(status)
        return status
    
    def check_sadi_consistence(self, fragment):
        '''
        check if same distance restraints make sense. Each length of an atom
        pair is tested agains the standard deviation of all distances.
        For a large standard deviation, the list is tested for outliers.
        :param atoms: atoms list of thr fragment
        :param restraints: restraints list
        :param fragment: frag name
        :param factor: factor for confidence interval
        '''
        atoms = self.db_dict[fragment]['atoms']
        restr = self.db_dict[fragment]['head']
        restraints = deepcopy(restr)
        atnames = [i[0].upper() for i in atoms]
        for num, line in enumerate(restraints):
            prefixes = []
            dev = 0.02
            line=line.split()
            if not line:
                continue
            if line[0].upper() == 'SADI':
                prefixes.append(line[0])
                del line[0]
                try:
                    if not str(line[0][0]).isalpha():
                        prefixes.append(line[0])
                        dev = line[0]
                        del line[0] # delete standard deviation
                except(IndexError):
                    return
                if len(line)%2 == 1: # test for uneven atoms count
                    print('Inconsistent SADI restraint line {} of "{}". Not all atoms form a pair.'.format(num, fragment))   
                pairs = pairwise(line)
                distances = []
                pairlist = []
                if len(pairs) <= 2:
                    return True
                for i in pairs:  
                    pairlist.append(i)
                    try:
                        a = atoms[atnames.index(i[0])][2:5]
                        b = atoms[atnames.index(i[1])][2:5]
                    except(ValueError):
                        return
                    a = [float(x) for x in a]
                    b = [float(y) for y in b]
                    dist = atomic_distance(a, b, self.get_unit_cell(fragment))
                    distances.append(dist)
                dist_minus_longest = sorted(distances)
                if len(distances) > 2:
                    del dist_minus_longest[-1]
                if len(distances) > 5:
                    del dist_minus_longest[-1]
                stdev = std_dev(dist_minus_longest)
                #print(stdev, fragment, stdev*len(pairs))
                # only do outlier test if standard deviation is suspiciously large:
                if stdev > 0.048:
                    outliers = nalimov_test(distances)
                    if outliers:
                        print("\nFragment {}:".format(fragment))
                        for x in outliers:
                            pair = ' '.join(pairlist[x])
                            print('Suspicious deviation in atom pair "{}" ({:4.3f} A, median: {:4.3f}) of SADI line {}:'.format(pair, distances[x], median(distances), num+1))
                            print(restr[num][:60], '...')
                            return False
                if stdev > 3*float(dev):
                    print("\nFragment {}:".format(fragment))
                    print('Suspicious restraints in SADI line {} with high Stdeviation {:4.3f} (median length: {:4.3f} A).'.format(num+1, stdev, median(distances)))
                    print(' '.join(prefixes+line))
        return True
                        
    
    def get_head_lines(self, fragment, db, line):
        '''
        return the head of the dbentry , the FRAG line and the comment of the 
        fragment as list of strings
        [['RESI CBZ', 'SADI C1 C2 C2 C3 C3 C4 C4 C5 C5 C6 C6 C1',
        'SADI Cl1 C2 Cl1 C6',
        'FLAT Cl1 > C6', 'SIMU Cl1 > C6', 'RIGU Cl1 > C6'], [FRAG 17 1 1 1 90 90 90], 
        [['Src:', 'pbe1pbe/6-311++G(3df,3pd),', 'Ilia', 'A.', 'Guzei'], 
        ['Name:', '1,2-Dichlorobenzene,', 'C6H4Cl2']] ]
        #####
        [list, list, dict]
        :param fragment: fragment name
        :type fragment:  string
        :param db:  database name e.g. dsr_db or dsr_user_db
        :type db: string
        :param line: line number where dbentry is located in the db file
        :type line: string
        :return nhead, # new head with unwrapped lines
             fragline, # line with frag command e.g. 'FRAG 17 1 1 1 90 90 90'
              comment: # comment lines from the head
        '''
        head = []
        nhead = []
        comment = []
        for i in self._db_plain_dict[db][int(line):]:
            if i.upper().startswith('FRAG'):  # collect the fragline
                fragline = i.rstrip(' \n\r')
                break
            head.append(i)
        for line in head:
            line = line.strip(' \n\r')
            line = ' '.join(line.split())
            if line.upper().startswith('REM'):
                comment.append(line.split()[1:])
                line = ''
                continue
            line = line.upper()
            nhead.append(line)
        # nhead is list of strings
        nhead = unwrap_head_lines(nhead)
        if not comment:
            comment = ['']
        try:
            if fragline:
                pass
        except(UnboundLocalError):
            print('Error. No cell parameters found in the database entry '\
                    'of "{}".'.format(fragment))
            print('Please add these parameters!')
            sys.exit(False)
        return (nhead, fragline, comment)
    
    def search_for_error_response(self, fragment):
        '''
        searches for a fragment name in the db as response to an invalid fragment name.
        :param fragment: the fragment
        :type fragment: string
        '''
        result = search_fragment_name(fragment)
        print('Do you mean one of these?:')
        print_search_results(result)

    def get_atoms_from_fragment(self, fragment):
        '''
        returns the atoms from the dbentry:
        [['O1', '1', '0.01453', '-1.6659', '-0.10966'],
        ['C1', '1', '0.00146', '-0.26814', '-0.06351'], ... ]
        :param fragment: fragment name
        :type fragment: string
        '''
        try:
            return self._dbentry_dict[fragment.lower()]['atoms']
        except KeyError:
            print('Could not find {} in database.'.format(fragment))
            self.search_for_error_response(fragment)
            sys.exit()

    def get_fragline_from_fragment(self, fragment):
        '''
        returns the line with FRAG 17 cell from the dbentry
        '''
        try:
            fragline = self._dbentry_dict[fragment.lower()]['fragline']
        except(KeyError):
            print('Fragment "{}" not found in database!'.format(fragment))
            self.search_for_error_response(fragment)
        return fragline

    def get_unit_cell(self, fragment):
        '''
        returns the unit cell parameters of the fragment
        [a, b, c, alpha, beta, gamma]
        '''
        return self.get_fragline_from_fragment(fragment)[2:]
    
    
    
    def get_line_number_from_fragment(self, fragment):
        '''
        returns the line number from the dbentry
        '''
        return self._dbentry_dict[fragment.lower()]['line']


    def get_head_from_fragment(self, fragment):
        '''
        returns the header of the dbentry of fragment.
        This header does not include comments, only the restraints.
        '''
        fragment = fragment.lower()
        try:
            head = self._dbentry_dict[fragment]['head']
        except KeyError:
            print('Could not find {} in database.'.format(fragment))
            self.search_for_error_response(fragment)
        self.check_db_header_consistency(fragment)
        self.check_sadi_consistence(fragment)
        return head


    def get_resi_from_fragment(self, fragment):
        '''
        returns the residue name of the dbentry of fragment
        can be either class or class + number.
        convention is only class.
        '''
        return self._dbentry_dict[fragment.lower()]['resi']


    def get_comment_from_fragment(self, fragment):
        '''
        returns the first comment line of the dbentry of a fragment
        if a line with "rem Name:" is present, this line is used as comment.
        :param fragment: actual fragment name
        :type fragment: string
        '''
        comment = self._dbentry_dict[fragment.lower()]['comment']
        for i in comment:
            if re.match(r'.*[n|N]ame:.*', i):
                i = i.split(' ', 1)[1:]
                comment = i
                break
        comment = ', '.join(comment)
        return comment

    def get_src_from_fragment(self, fragment):
        '''
        returns the source line of the dbentry of a fragment
        if a line with "rem Src:" is present.
        :param fragment: actual fragment name
        :type fragment: string
        '''
        src = self._dbentry_dict[fragment.lower()]['comment']
        for i in src:
            if re.match(r'.*[s|S][r|R][c|C]|Source:.*', i):
                i = i.split(' ', 1)[1:]
                src = i
                break
        src = ', '.join(src)
        return src


    def get_db_name_from_fragment(self, fragment):
        '''
        returns the fragment database name of fragment x
        '''
        return self._dbentry_dict[fragment.lower()]['db']


class ImportGRADE():

    def __init__(self, grade_tar_file, invert=False):
        '''
        class to import fragments from GRADE of Global Phasing Ltd.
        :param gradefile:
        :type gradefile:
        :param invert:
        :type invert:
        '''
        try:
            main_dbdir = os.environ["DSR_DIR"]
        except(KeyError):
            main_dbdir='./'
        self.el = Element()
        self.invert = invert
        self._getdb = ReadDB(main_dbdir)
        self._db_dir = expanduser("~")
        self._db_tags = self._getdb.find_db_tags()
        self._gdb = global_DB(invert=False)
        self._db = self._gdb.build_db_dict()
        gradefiles = self.get_gradefiles(grade_tar_file)
        self._pdbfile = gradefiles[0]
        self._dfixfile = gradefiles[1]
        #self._obpropfile = gradefiles[2]
        self._atoms = self.get_pdbatoms(self._pdbfile)
        self._firstlast = self.get_first_last_atom(self._atoms)
        self._restraints = self.get_restraints()
        self._resi_name = self.get_name_from_pdbfile(self._pdbfile)
        if not isinstance(self._resi_name, str):
                self._resi_name = self._resi_name.decode()
        self._comments = self.get_comments()


    def get_gradefiles(self, grade_tar_file):
        '''
        returns the .mol2 and .dfix file location
        '''
        grade_base_filename = os.path.splitext(grade_tar_file)
        if grade_base_filename[1] == '.tgz':
            try:
                gradefile = tarfile.open(grade_tar_file)#, encoding="ascii")
            except(IOError):
                print('No such file or directory: {}'.format(grade_tar_file))
                sys.exit(0)
        else:
            print('File {} is not a valid file to import from!'.format(grade_base_filename[1]))
            sys.exit(0)
        pdbfile = False
        dfixfile = False
        for i in gradefile.getnames():
            if i.endswith('.pdb'):
                if re.match(r'.*obabel.*', i):
                    continue
                pdbfile = gradefile.extractfile(i)
            if i.endswith('.dfix'):
                dfixfile = gradefile.extractfile(i)
        if not pdbfile:
            print(' .pdb file not found!')
            self.import_error(grade_tar_file)
        elif not dfixfile:
            print(' .dfix file not found!')
            self.import_error(grade_tar_file)
        output = []
        output.append(self._read_tarfile_as_ascii(pdbfile))
        output.append(self._read_tarfile_as_ascii(dfixfile))
        return output


    def _read_tarfile_as_ascii(self, filehandle):
        tmp = []
        for line in filehandle:
            if not isinstance(line, str):
                line = line.decode('ascii')
            tmp.append(line)
        return tmp


    def get_name_from_pdbfile(self, pdbfile):
        '''
        get the fragment name from the pdbfile.txt file
        :param pdbfile: file with some information about the molecule
        :type pdbfile: list of strings
        '''
        regex = re.compile(r'^.*Compound full name.*')
        for line in pdbfile:
            if not isinstance(line, str):
                line = line.decode('ascii')
            if regex.match(line):
                line = line.replace('_', '')
                line = line.replace('-', '')
                line = line.replace('#', '')
                break
        if not line:
            return 'NONE'
        try:
            line = line.split()
            line[4]
        except(IndexError):
            return 'NONE'
        return line[5]


    def get_comments(self):
        '''
        returns a detailed comment out of the .dfix file where the fragment was imported from

        #   REM Produced by Grade Web Server http://grade.globalphasing.org
        #   REM GEN: Generated by GRADE 1.2.5 (December 20 2013)
        #   REM GEN: from SMILES C1=CC=C2C(=C1)C=C(C3=CC=CC=C23)Br
        #   REM GEN: using MOGUL 1.6(RC5), CSD as535be with quantum mechanics RM1
        #   REM grade-cif2shelx output
        #   REM Version: 0.0.5 <Dec 20 2013>
        '''
        matches = ['REM Produced by Grade', 'REM GEN:', 'REM grade-cif2shelx', 'REM Version:', 'REM Total charge']
        comments = []
        name = 'REM Name: '+self._resi_name
        comments.append(name.split())
        for m in matches:
            for line in self._dfixfile:
                if not isinstance(line, str):
                    line = line.decode()
                if line.startswith(m):
                    comments.append(line.split())
        return comments


    def get_first_last_atom(self, atoms):
        '''
        returns the first and the last atom from the imported atom list
        '''
        try:
            first = atoms[0][0]
            last = atoms[-1][0]
        except:
            return False
        return (first, last)


    def get_restraints(self):
        '''
        reads the restraints from a dfix-file
        '''
        restraints = []
        for line in self._dfixfile:
            if not isinstance(line, str):
                line = line.decode()
            line = line.strip('\n\r').split()
            if line:
                line[0] = line[0][:4].upper()
                if line[0] in RESTRAINT_CARDS:
                    restraints.append(line)
        if self._firstlast:
            atoms = ' > '.join(self._firstlast)
            atoms_s = 'SIMU ' + atoms
            atoms_s = atoms_s.split()
            atoms_r = 'RIGU ' + atoms
            atoms_r = atoms_r.split()
            restraints.append(atoms_s)
            restraints.append(atoms_r)
        return restraints


    def get_pdbatoms(self, pdb):
        '''
        returns the atoms from a pdb file
        '''
        atomlines = []
        for line in pdb:
            if not isinstance(line, str):
                line = line.decode()
            line = line.split()
            if line[0] == 'HETATM':
                if line[-1] == 'H':
                    continue
                tmp = []
                tmp.append(line[2])
                tmp.append(line[11])
                tmp.extend(line[6:9])
                atomlines.append(tmp)
        if not atomlines:
            self.import_error(pdb)
        if self.invert:
            atomlines = invert_dbatoms_coordinates(atomlines)
        return atomlines


    def bild_grade_db_entry(self):
        '''
        builds a dbentry from the information supplied by GRADEs
        .mol2 and .dfix file
        '''
        db_import_dict = {}
        num = 1
        name = self._resi_name[:3].upper()
        if not isinstance(name, str):
            name = name.decode()
        resi_name =  name + str(num)
        if not self._db_tags:
            print('Unable to import fragment. Database is empty.')
            sys.exit(False)
        for i in self._db_tags:
            while resi_name.upper() == i[0]:
                num = num + 1
                resi_name = resi_name[:3] + str(num)
        # print 'using {} as resiname'.format(resi_name)
        fragline = 'FRAG 17 1  1  1  90  90  90'
        db_import_dict[resi_name] = {
                'head'    : self._restraints,
                'resi'    : resi_name,
                'fragline': fragline.split(),
                'atoms'   : self._atoms,
                'line'    : None,
                'db'      : 'dsr-user-db',
                'comment' : self._comments,
                'name'    : resi_name
                }
        return db_import_dict


    def import_error(self, filename):
        '''
        warns for import errors
        '''
        print('Unable to import GRADE file {}'.format(filename))
        print('GRADE import relies on GRADE v1.100 and up.')
        sys.exit(False)


    def write_user_database(self):
        '''
        writes content of existing dsr_user_db.txt and the imported GRADE entry to
        the dsr_user_db.txt
        '''
        filename = os.path.join(self._db_dir, 'dsr_user_db.txt')
        fragments = list(self._db.keys())
        imported_entry = self.bild_grade_db_entry()
        grade_db_names = list(imported_entry.keys())
        user_db_names = []
        for i in fragments:
            if self._db[i]['db'] == 'dsr_user_db':
                user_db_names.append(self._db[i]['name'])
        # try to write the new dbentry:
        try:
            with open(filename, 'w') as f:
                for name in grade_db_names:
                    print('Importing {} ({}) to user database...'.format(self._resi_name, name))
                    atomlist = imported_entry[name]['atoms']
                    comment = imported_entry[name]['comment']
                    comment = '\n'.join([' '.join(i) for i in comment if i])
                    head = '\n'.join([' '.join(x) for x in imported_entry[name]['head']])
                    atoms = '\n'.join(['{:<6} -{}  {:>8.3f}{:>8.3f}{:>8.3f}'\
                            .format(y[0], self.el.get_atomic_number(y[1]),float(y[2]), float(y[3]), float(y[4])) for y in atomlist])
                    resi_name = str(name)
                    cell = '  '.join(imported_entry[name]['fragline'])
                    dbentry = '<{}> \n{} \nRESI {} \n{} \n{} \n{} \n</{}>\n''\
                        '.format(resi_name, comment, resi_name, head, cell, atoms, resi_name)
                    f.write(dbentry)
        except(IOError) as e:
            print(e)
            sys.exit(-1)
        # try to write existing dbentries:
        try:
            with open(filename, 'a+') as fu:
                for i in fragments:
                    name = i
                    if self._db[i]['db'] == 'dsr_user_db':
                        # userdb = list(self._db[i].keys())
                        atomlist = self._db[i]['atoms']
                        head = '\n'.join([''.join(x) for x in self._db[i]['head']])
                        atoms = '\n'.join(['{:<6}{:<2}{:>8.3f}{:>8.3f}{:>8.3f}'\
                            .format(y[0], y[1], float(y[2]), float(y[3]), float(y[4])) for y in atomlist])
                        resi_name = self._db[i]['resi']
                        comment = self._db[i]['comment']
                        comment = '\nREM '.join(comment)
                        fragline = '  '.join(self._db[i]['fragline'])
                        dbentry = '\n<{}> \nREM {} \nRESI {} \n{} \n{} \n{} \n</{}>\n'\
                            ''.format(name, comment, resi_name, head, fragline, atoms, name)
                        fu.write(dbentry)
        except(IOError) as e:
            print(e)
            sys.exit(-1)
        print('User database successfully updated.')




if __name__ == '__main__':

    gdb = global_DB(invert=False)
    db = gdb.build_db_dict()
    dbnames = list(db.keys())
    
    
    for names in dbnames:
        # fragment = 'pfanion'
        fragment = names    
        # fragline = gl.get_fragline_from_fragment(fragment)  # full string of FRAG line
        # dbatoms = gl.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
        dbhead = gdb.get_head_from_fragment(fragment)  # this is only executed once
        dbhead = unwrap_head_lines(dbhead)
        dbatoms = gdb.get_atoms_from_fragment(fragment)
        # print dbatoms
        # print('residue:', db['toluene']['resi'])
        # print('line of db:', db['toluene']['line'])
        # print('database:', db['toluene']['db'])
        # print(fragline)
    
        # for i in dbatoms:
        #    print i
        head = db[fragment]['head']
        gdb.check_sadi_consistence(fragment)
    sys.exit()
    # mog = ImportGRADE('./test-data/ALA.gradeserver_all.tgz')
    # mog = ImportGRADE('./test-data/LIG.gradeserver_all.tgz')

    # import tempfile
    # import tarfile
    gf = './test-data/LIG.gradeserver_all.tgz'
    gradefile = tarfile.open(gf)
    # create a temporary file and download platon to it.

#    print
    # list content of tgz file
    # print gradefile.getmembers()

    for i in gradefile.getnames():
        if i.endswith(('.pdb', '.dfix', '.mol2')):
            if re.match('.*obabel.*', i):
                continue
            # localFile = tempfile.TemporaryFile()
            print(i)
            gfile = gradefile.extractfile(i)
            # localFile.write(gfile)
#            print gfile.readlines()
            # print i



    # print localFile.readlines()
    print()
    # print misc.ll_to_string(mog.get_molatoms())
    # for i in mog.get_molatoms('./test-data/TOL.mol2'):
    #    #print '{:<6}{:>1}{:>10}{:>10}{:>10}'.format(i[0], i[1], i[2], i[3], i[4])
    #    print '{:<6}{:>1}{:>10}{:>10}{:>10}'.format(*i)
    #    #print i[3], i[0], i[1], i[2]
    # print mog.get_atomnumbers()

    # for i in mog.get_restraints('./test-data/TOL.dfix'):
    #    print i

    # mog.bild_grade_db_entry()
    # print mog.bild_grade_db_entry('{}.mol2', '{}.dfix').format(dsr_dict[import_grade], dsr_dict[import_grade])
    # mog.write_user_database()
