# -*- encoding: utf-8 -*-
# möp
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <dkratzert@gmx.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#


import os
import re
import sys
import tarfile
from copy import deepcopy

from src.atomhandling import get_atomtypes
from src.atoms import Element
from src.constants import atomregex, SHX_CARDS, RESTRAINT_CARDS, sep_line
from src.misc import atomic_distance, nalimov_test, std_dev, median, pairwise, \
    unwrap_head_lines, dice_coefficient2, frac_to_cart

not_existing_error = '*** Fragment "{}" not found in database ***'


def invert_atomic_coordinates(atoms):
    """
    Inverts SHELXL atom coordinates.

    :param atoms: list of atom list

    >>> atoms = [['C1', 1, 0.44, 0.21, -1.23 ], ['C2', 1, -1.44, -0.21, 21.23]]
    >>> invert_atomic_coordinates(atoms)
    [['C1', 1, -0.44, -0.21, 1.23], ['C2', 1, 1.44, 0.21, -21.23]]
    >>> invert_atomic_coordinates([['C2', 1, 'd', 1, 2]])
    Unable to invert fragment coordinates.
    []
    """
    catoms = deepcopy(atoms)
    for line in catoms:
        try:
            inv_coord = [-x for x in line[-3:]]
        except:
            print('Unable to invert fragment coordinates.')
            return []
        line[-3:] = inv_coord
    return catoms


def search_fragment_name(search_string, gdb, numresults=6):
    """
    searches the Name: comments in the database for a given name

    :param numresults: number of results to return after search
    :type search_string: str
    :param gdb: database object
    :type gdb: global_DB
    """
    names_list = []
    for fragment in gdb.databases:
        fragname = gdb[fragment]['name']
        line_number = gdb[fragment]['startline']
        dbname = gdb[fragment]['dbname']
        names_list.append([fragment, fragname, line_number, dbname])
    search_results = []
    for fragment in names_list:
        key = make_sortkey(fragment[1], searchkey=True)
        # key[0] and key[1] combined to have also the sum formula:
        fullstring  = "{0}{1}".format(key[0], key[1])
        if fragment[1].lower().startswith(search_string):
            # give names starting with a certain string high priority
            coefficient = 3.0
        elif search_string in fullstring or search_string in fragment[0]:
            coefficient = 2.1
        else:
            coefficient = dice_coefficient2(search_string, fullstring)
        fragment.append([coefficient, key[1]])
        if coefficient > 0.1:
            search_results.append(fragment)
    # select the best n results, sort for first and second search key
    selected_results = sorted(search_results, key=lambda coeff: [coeff[-1][0], coeff[-1][1]], reverse=True)[
                       :numresults]
    return selected_results


def print_search_results(results):
    """
    prints the results of a database search to screen and exit.
    results are
    """
    print(' Fragment          | Full name, Comments                      | Line number')
    print(sep_line)
    for line in results:
        print(' {:15s}   | {:40s} | {}'.format(line[0], line[1], line[2]))
    sys.exit()


def make_sortkey(full_name, searchkey=False):
    """
    Algorythm inspired by W. Sage J. Chem. Inf: Comput. Sci. 1983, 23, 186-197
    """
    keylist = []
    full_name = ''.join(e for e in full_name if e not in ('{}()[],'))
    full_name = full_name.split(' ')
    rest = " " + ' '.join(full_name[1:])
    full_name = full_name[0].lower()
    numbers = ''.join(e for e in full_name if e in r'0123456789')
    if full_name.startswith('tert-'):
        full_name = full_name[4:]
    if full_name.startswith('sec-'):
        full_name = full_name[3:]
    if full_name.startswith('iso-'):
        full_name = full_name[3:]
    if full_name.startswith('bis-'):
        full_name = full_name[3:]
    if full_name.startswith('tris-'):
        full_name = full_name[3:]
    if full_name.startswith('mono-'):
        full_name = full_name[4:]
    if full_name.startswith('t-'):
        full_name = full_name[1:]
    if full_name.startswith('p-'):
        full_name = full_name[2:]
    if full_name.startswith('o-'):
        full_name = full_name[1:]
    if full_name.startswith('n-'):
        full_name = full_name[1:]
    if full_name.startswith('nn-'):
        full_name = full_name[1:]
    if full_name.startswith('nnn-'):
        full_name = full_name[1:]
    if full_name.startswith('m-'):
        full_name = full_name[1:]
    if full_name.startswith('i-'):
        full_name = full_name[1:]
    full_name = ''.join(e for e in full_name if e not in ('+-_.\'1234567890, '))
    if searchkey:
        keylist = [full_name + rest, numbers]  # enables search for sum formulae
    else:
        keylist = [full_name, numbers]
    return keylist





# hardwired names of the database files:
# dsr_db.txt is the db from the distribution. This file should not be edited.
# dsr_user_db.txt if this file exists, all its content is also read in.

def read_file_data(filepath):
    # type: (str) -> list
    """
    reads the database files and returns them as list.
    """
    dblist = []
    try:
        with open(filepath, 'r') as f:
            dblist = f.read().splitlines()
    except (IOError, TypeError) as e:
        print(e)
    return dblist


class ParseDB(object):
    """
    reads in the system and user db files and makes
    a dictionary of them.
    maindb = read_file_data(path_to_maindb)
    """

    def __init__(self, maindb_path=None, userdb_path=None):
        # type: (str, str) -> NotImplemented
        """
        self._databases: dictionary with the individial fragments
        """
        self.databases = {}
        self.maindb_path = maindb_path
        self.userdb_path = userdb_path
        if maindb_path:
            self.parse(maindb_path, 'dsr_db')
        if userdb_path:
            self.parse(userdb_path, 'dsr_user_db')

    def parse(self, dbpath='./dsr_db.txt', dbname='dsr_db'):
        # type: (str, str) -> dict
        """
        This method returns all fragment name tags in the database

        >>> dbpath = os.path.abspath('./dsr_db.txt')
        >>> db = ParseDB(dbpath)
        >>> db.parse(dbpath, 'dsr_db')['water']['name'] # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        'Water, H2O'
        >>> db.parse(dbpath, 'dsr_db')['water']['comments']
        []
        >>> db.parse(dbpath, 'dsr_db')['water']['atoms'] # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        [['O1', 4, 0.0, 0.0, 0.0],
        ['H1', 2, 0.9584, 0.0, 0.0],
        ['H2', 2, -0.2392, 0.9281, 0.0]]
        >>> db.parse(dbpath, 'dsr_db')['water']['cell'] # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        [1.0, 1.0, 1.0, 90.0, 90.0, 90.0]
        >>> db.parse(dbpath, 'dsr_db')['water']['source'] # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        'pbe1pbe/6-311++G(3df,3pd), Ilia A. Guzei'
        >>> db.parse(dbpath, 'dsr_db')['water']['resi'] # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        'H2O'
        >>> db.parse(dbpath, 'dsr_db')['water']['restraints'] # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        ['DFIX 0.9584 0.001 O1 H1 O1 H2', 'DFIX 1.5150 0.001 H1 H2']
        >>> db.parse(dbpath, 'dsr_db')['water']['dbname'] # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        'dsr_db'
        >>> db.parse(dbpath, 'dsr_db')['water']['hfix'] # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        []
        """
        frag_tag = ''
        db = {}
        start_regex = re.compile(r'<[^/].*>', re.IGNORECASE)  # regular expression for db tag.
        starttag = False
        fraglines = []
        end_regex = None
        startnum = 0
        for num, line in enumerate(read_file_data(dbpath)):
            if line.startswith('#'):
                continue
            if not line:
                continue
            # matching end tag
            if end_regex and end_regex.match(line):
                starttag = False
                db[frag_tag].update(
                        {'endline'  : num + 1,
                         'startline': startnum + 1})
                db = self.parse_fraglines(frag_tag, fraglines, db)
                fraglines = []
            # start tag was found, appending lines to fragment list
            if starttag:
                fraglines.append(line)
            # matching start tag and compiling end regex
            if start_regex.match(line):
                if starttag:
                    print('*** Error in database "{}.txt" in line {}. End tag is missing ***'.format(dbname, num + 1))
                # lower case is essential here:
                frag_tag = line.strip('<> \n\r').lower()
                if frag_tag in db:
                    print('\n*** Duplicate database entry "{}" found! Please remove/rename '
                          'second entry\nand/or check all end tags in the database dsr_usr_db.txt '
                          'or dsr_db.txt. ***'.format(frag_tag))
                    sys.exit()
                starttag = True
                startnum = num
                db[frag_tag] = {'dbname': dbname}
                end_regex = re.compile(re.escape(r'</{}>'.format(frag_tag)), re.IGNORECASE)
                continue
        self.databases.update(db)
        return db

    @staticmethod
    def parse_fraglines(fragname_tag, fraglines, db):
        # type: (str, list, dict) -> dict
        """
        Fills the database dictionary with fragment data.
        """
        fragname_tag = fragname_tag.lower()
        headlist = []
        comments = []
        hfix = []
        residue = ''
        atoms = []
        cell = []
        source = ''
        name = ''
        nameregex = re.compile(r'REM\s+NAME:', re.IGNORECASE)
        srcregex = re.compile(r'REM\s+(SRC:|SOURCE:)', re.IGNORECASE)
        # devide atoms and the rest:
        for num, aline in enumerate(fraglines):
            if atomregex.match(aline):  # search atoms
                atline = aline.split()[:5]  # convert to list and use only first 5 columns
                if atline[0] not in SHX_CARDS:  # exclude all non-atom cards
                    coords = []
                    try:
                        coords = [float(x) for x in atline[2:]]
                        atline[1] = int(atline[1])
                    except ValueError:
                        print("*** Invalid atomic coordinates in line {} of {}.txt (Fragment: {}) ***"
                              .format(db[fragname_tag]['startline'] + num + 1, db[fragname_tag]['dbname'],
                                      fragname_tag))
                        sys.exit()
                    atline[2:] = coords
                    atoms.append(atline)
                    fraglines[num] = ''
        # bring all wrapped lines with = at end to a single line:
        fraglines = unwrap_head_lines(fraglines)
        for num, line in enumerate(fraglines):
            # collect the comments:
            if line.upper().startswith('REM'):
                if line.upper().startswith('REM HFIX'):
                    hfix.append(line)
                    continue
                # Where do the fragment come from?
                if srcregex.match(line):
                    source = ' '.join(line.split()[2:])
                    fraglines[num] = ''
                    continue
                # Get the full name:
                if nameregex.match(line):
                    name = ' '.join(line.split()[2:])
                    fraglines[num] = ''
                    continue
                comments.append(line)
                continue
            command = line[:4].upper().strip()
            if command and command in SHX_CARDS:
                if command == 'RESI':
                    resiline = line.split()
                    for n, i in enumerate(resiline):
                        if not i[0].isalpha():
                            del resiline[n]
                    try:
                        residue = resiline[1].upper()
                        continue
                    except IndexError:
                        print('*** Invalid residue definition in database entry {}. ***'.format(fragname_tag))
                        sys.exit()
                # collect the unit cell of the fragment:
                if command == 'FRAG':
                    try:
                        cell = [float(x) for x in line.split()[2:8]]
                    except (ValueError, IndexError):
                        print('*** Invalid unit cell found in database entry {}. ***'.format(fragname_tag))
                        sys.exit()
                    continue
                # these must be restraints:
                if line[:4] in RESTRAINT_CARDS:
                    headlist.append(line)
            elif command:
                print('*** Bad line {} in database entry "{}" found! ({}.txt) ***'
                      .format(num + db[fragname_tag]['startline'] + 1, fragname_tag, db[fragname_tag]['dbname']))
                print(line)
        if not cell:
            print('*** Error. No cell parameters or malformed cell found in the database entry ' \
                  'of "{}" ***'.format(fragname_tag))
            sys.exit()
        if not name:
            # This happens if no name is given or "REM Name:" has errors (The name is explicitely not mandatory!).
            name = fragname_tag
        if not db:
            print('*** No database found! ***\n')
            sys.exit()
        if not atoms:
            # Can not print this out, because shelXle GUI will fail if text is printed.
            # print('*** No atoms found in database entry {} line {} of {}.txt***'.format(fragname_tag,
            #                                    db[fragname_tag]['startline'] + 1, db[fragname_tag]['dbname']))
            del db[fragname_tag]
            return db
        db[fragname_tag].update({
            'restraints': headlist,  # header with just the restraints
            'resi'      : residue,  # the residue class
            'cell'      : cell,  # FRAG ...
            'atoms'     : atoms,  # the atoms as lists of list
            'comments'  : comments,  # the comment line
            'source'    : source,
            'name'      : name,
            'hfix'      : hfix,
        })
        return db

    def __getitem__(self, fragment):
        # type: (str) -> dict
        try:
            return self.databases[fragment]
        except KeyError:
            print("*** Fragment {} was not found in the Database! ***".format(fragment))
            sys.exit()

    def __iter__(self):
        """
        >>> db = ParseDB('./dsr_db.txt')
        >>> [x for x in db][:3]
        ['12-dichlorobenz', '12-difluorobenz', '12c4']
        """
        # return iter(self.databases.keys()) #python2
        return iter(list(sorted(self.databases.keys())))

    def list_fragments(self):
        # type: () -> list
        """
        list all fragments in the db as list of lists
        [['tbu-c', 1723, 'dsr_db', 'Tert-butyl-C'], ...]
        """
        fraglist = []
        for frag in self.databases:
            name = self.get_fragment_name(frag)
            key = make_sortkey(name)
            line = [frag,
                    self.get_startline(frag),
                    self.get_db_name(frag),
                    name,
                    key[0],
                    key[1]]
            fraglist.append(line)
        fraglist.sort(key=lambda x: (x[4], x[5]))
        return fraglist

    def get_atomic_numbers(self, fragment=''):
        # type: (str) -> list
        """
        Returns the atomic numbers of the atoms in a fragment in
        same order as in the database as list.
        """
        atnumbers = []
        types = get_atomtypes(self.get_atoms(fragment))
        el = Element()
        for i in types:
            atnumbers.append(el.get_atomic_number(i))
        return atnumbers

    def get_atomnames(self, fragment='', uppercase=False):
        # type: (str, bool) -> list
        """
        Returns a list of only the atom names of a fragment.
        """
        atoms = self.get_atoms(fragment)
        names = []
        for i in atoms:
            if uppercase:
                names.append(i[0].upper())
            else:
                names.append(i[0])
        return names

    def get_head_for_gui(self, fragment):
        """
        returns header information of the specific fragment:
        tag, Name/comment, source, cell, residue, dbtype, restr, atoms
        >>> dbpath = os.path.abspath('./dsr_db.txt')
        >>> db = ParseDB(dbpath)
        >>> db.get_head_for_gui('benZene')
        ... # doctest: +NORMALIZE_WHITESPACE
        <tag>
         benzene
        </tag>
        <comment>
         Benzene, Phenyl, C6H6
        </comment>
        <source>
         UGEDEQ
        </source>
        <cell>
         1;;1;;1;;90;;90;;90
        </cell>
        <residue>
         BENZ
        </residue>
        <dbtype>
         dsr_db
        </dbtype>
        <restr>
         SADI 0.02 C1 C2 C2 C3 C3 C4 C4 C5 C5 C6 C6 C1;;SADI 0.04 C1 C5 C1 C3 C5 C3 C4 C2 C4 C6 C6 C2;;FLAT C1 > C6;;SIMU C1 > C6;;RIGU C1 > C6
        </restr>

        """
        fragment = fragment.lower()
        print("<tag>\n", fragment, "\n</tag>")
        print("<comment>\n", self.get_fragment_name(fragment), "\n</comment>")
        print("<source>\n", self.get_src(fragment), "\n</source>")
        # it returns cartesian coordinates, so we need this cell here:
        print("<cell>\n", '1;;1;;1;;90;;90;;90', "\n</cell>")
        print("<residue>\n", self.get_resi(fragment), "\n</residue>")
        print("<dbtype>\n", self.get_db_name(fragment), "\n</dbtype>")
        print("<restr>\n", ';;'.join(self.get_restraints(fragment)), '\n</restr>')
        self.check_consistency(fragment)
        self.check_db_restraints_consistency(fragment)
        if not self.check_sadi_consistence(fragment):
            sys.exit()
        self.check_db_atom_consistency(fragment)

    def check_consistency(self, fragment):
        # type: (str) -> bool
        """
        check if the fragline makes sense and if the fragment_dict is complete
        """
        fragment = fragment.lower()
        if fragment in ['cf3', 'cf6', 'cf9']:
            return True
        try:
            dbentry = self.databases[fragment]
        except KeyError:
            print("*** Fragment {} not found in the database. ***".format(fragment))
            sys.exit()
        if not dbentry['cell']:
            print('*** Error. No cell parameters or malformed cell found in the database entry ' \
                  'of "{}" ***'.format(fragment))
            print('*** Please add these parameters! ***')
            sys.exit()
        if not dbentry['restraints']:
            print('*** Restraints in the header of database entry "{}" ({}) missing! Check your Database! ***' \
                  .format(fragment, dbentry['name']))
            return False
        if not dbentry['resi']:
            print('*** Residue in the header of database entry "{}" ({}) missing! Check your Database! ***' \
                  .format(fragment, dbentry['name']))
        if not dbentry['atoms']:
            print('*** No atoms found in database entry {} line {} of {}.txt***'.format(fragment,
                                                                                        dbentry['startline'] + 1,
                                                                                        dbentry['dbname']))
            print('*** Have you really followed the syntax? ***')
            sys.exit()
        if not dbentry['endline']:
            print('*** Could not find end of dbentry for fragment "{}" in line {} of "{}.txt". '
                  'Check your database. ***'.format(fragment, dbentry['startline'] + 1, dbentry['dbname']))
            sys.exit()
        return True

    def get_sum_formula(self, fragment):
        """
        returns the sum formula of a fragment
        as string like 'SI4 C9 O1'
        :param fragment: fragment name
        :type fragment: string
        """
        fragment = fragment.lower()
        try:
            types = get_atomtypes(self.databases[fragment]['atoms'])
        except:
            return None
        formula = ''
        for el in set(types):
            num = types.count(el)
            formula += '{}{} '.format(el, num)
        return formula

    def check_db_atom_consistency(self, fragment):
        # type: (str) -> bool
        """
        This method is for atoms only (without db header)!
        check the db for duplicates:
        """
        dbatoms = self.get_atomnames(fragment, uppercase=True)
        # check for duplicates:
        while dbatoms:
            at = dbatoms.pop()
            if at in dbatoms:
                print('*** Duplicate atom {0} in database entry "{1}" ({2}) '
                      'found! Check your database... ***'
                      .format(at, fragment, self.databases[fragment]['name']))
                sys.exit()
        return True

    def get_hfixes(self, fragment, resi_class):
        """
        Returns pre-defined hfix instructions.
        :param fragment: fragment name
        :return: A string with the REM HFIX ... instruction
        >>> from src.misc import wrap_headlines
        >>> hf = ['REM HFIX 123 C1 C2 C3 C4 C5']
        >>> wrap_headlines(hf, 18)
        ['REM HFIX 123 C1 C2 =\\n   C3 C4 C5\\n']
        >>> db = ParseDB('./dsr_db.txt')
        >>> db.get_hfixes('Adamantane', 'ADAM')
        ['REM HFIX_ADAM 23 C2 C3 C4 C6 C8 C10', 'REM HFIX_ADAM 13 C5 C7 C9']
        >>> db.get_hfixes('Adamantane', '')
        ['REM HFIX 23 C2 C3 C4 C6 C8 C10', 'REM HFIX 13 C5 C7 C9']
        >>> db.get_restraints('Adamantane')
        ['SADI C1 C2 C1 C3 C1 C4 C5 C6 C6 C7 C7 C8 C8 C9 C9 C10 C10 C5 C2 C9 C3 C7 C4 C5', 'SADI 0.04 C2 C3 C3 C4 C4 C2 C5 C9 C9 C7 C7 C5 C6 C10 C6 C8 C8 C10', 'SADI 0.05 C5 C8 C7 C10 C9 C6', 'SIMU C1 > C10', 'RIGU C1 > C10']
        >>> db.get_hfixes('oc(cf3)3', '')
        []
        """
        fragment = fragment.lower()
        try:
            hfix = [' '.join(x.upper().split()) for x in self.databases[fragment]['hfix']]
        except KeyError:
            print(not_existing_error.format(fragment))
            sys.exit()
        if hfix and resi_class:
            resihfix = []
            for line in hfix:
                # Adding _resiclass to HFIX instructions:
                # REM HFIX C1 C2  ->  REM HFIX_class C1 C2 
                line = '{}_{} {}'.format(line[:8], resi_class, line[9:])
                resihfix.append(line)
            return resihfix
        return hfix

    def check_db_restraints_consistency(self, fragment):
        """
        - Checks if the Atomnames in the restraints of the dbhead are also in
          the list of the atoms of the respective dbentry.
        - Checks wether restraints cards are vaid.
        """
        fragment = fragment.lower()
        atoms = self.get_atoms(fragment)
        status = True
        atoms = [i[0].upper().strip() for i in atoms]
        restraint_atoms_list = set([])
        for line in self.get_restraints(fragment):
            for i in line.split():
                # Todo: resolve ranges:
                if i in ('>', '<'):
                    continue
                try:
                    # Test if first parameter is a distance or a standard deviation:
                    if i in ['DFIX', 'DANG']:
                        # print('####')
                        if len(line.split()) > 1:  # there is more than just DFIX
                            # Test if this is correct: DFIX d s[0.02] atom pairs (d should be greater s)
                            if float(line.split()[1]) < float(line.split()[2]):
                                print('\n*** Sigma > d in DFIX/DANG of "{}: {}" '
                                      'in database {}. ***'
                                      .format(fragment, self.get_fragment_name(fragment), self.get_db_name(fragment)))
                        else:
                            print('\n*** Incomplete DFIX/DANG in restraints of "{}: {}" in database {}. ***'
                                  .format(fragment, self.get_fragment_name(fragment), self.get_db_name(fragment)))
                    # if i in SHX_CARDS:
                    #    continue
                    else:
                        # just test if there is not a number:
                        float(i)
                except ValueError:
                    if i in SHX_CARDS:
                        continue
                    restraint_atoms_list.add(i)
        for atom in restraint_atoms_list:
            if atom.upper() not in atoms:
                status = False
                print('\n*** Unknown atom "{}" in restraints of "{}: {}" in database {}. ***'
                      .format(atom, fragment, self.get_fragment_name(fragment), self.get_db_name(fragment)))
        if not status:
            print('*** Check database entry. ***\n')
            sys.exit()
        return status

    def check_sadi_consistence(self, fragment):
        """
        check if same distance restraints make sense. Each length of an atom
        pair is tested agains the standard deviation of all distances.
        For a large standard deviation, the list is tested for outliers.
        :param fragment: frag name
        """
        atoms = self.get_atoms(fragment)
        restr = self.get_restraints(fragment)
        restraints = deepcopy(restr)
        atnames = self.get_atomnames(fragment, uppercase=True)
        good = True
        for num, line in enumerate(restraints):
            prefixes = []
            dev = 0.02
            line = line.split()
            if not line:
                continue
            if line[0].upper() == 'SADI':
                prefixes.append(line[0])
                del line[0]
                try:
                    if not str(line[0][0]).isalpha():
                        prefixes.append(line[0])
                        dev = line[0]
                        del line[0]  # delete standard deviation
                except IndexError:
                    return False
                if len(line) % 2 == 1:  # test for uneven atoms count
                    print('*** Inconsistent SADI restraint line {} of "{}". '
                          'Not all atoms form a pair ***'.format(num, fragment))
                pairs = pairwise(line)
                distances = []
                pairlist = []
                if len(pairs) <= 1:
                    return True
                for i in pairs:
                    if i in pairlist or tuple(reversed(i)) in pairlist:
                        print('*** Duplicate atom pair "{}" in SADI restraint line {} of "{}". ***'.format(" ".join(i),
                                                                                                           num,
                                                                                                           fragment))
                    pairlist.append(i)
                    try:
                        atom1 = [float(x) for x in atoms[atnames.index(i[0])][2:5]]
                        atom2 = [float(x) for x in atoms[atnames.index(i[1])][2:5]]
                    except ValueError:
                        return False
                    dist = atomic_distance(atom1, atom2, self.get_cell(fragment))
                    distances.append(dist)
                if len(pairlist) == 2:
                    # Find restraints with one pair, where 1,2 and 1,3 distances are mixed:
                    pairdev = 1.0 - (min(distances) / max(distances))
                    if pairdev > (2.0 * float(dev)):
                        print('*** Suspicious deviation of {:.3f} A for "{}" in {} ***'.format(pairdev, restraints[num],
                                                                                               fragment))
                        return False
                    return True  # because esd is too sensitive for two distances.
                stdev = std_dev(distances)  # Error distribution of
                # only do outlier test if standard deviation is suspiciously large:
                if stdev > 0.065:
                    outliers = nalimov_test(distances)
                    if outliers:
                        print("\nFragment {}:".format(fragment))
                        for x in outliers:
                            pair = ' '.join(pairlist[x])
                            print('*** Suspicious deviation of atom pair "{}" ({:4.3f} A, median: {:4.3f}) after '
                                  'line {} in {}.txt ***'.format(pair, distances[x], median(distances),
                                                                 self.get_startline(fragment) + num + 1,
                                                                 self.get_db_name(fragment))
                                  )
                            print('*** {} ... ***'.format(restr[num][:60]))
                            good = False
                if (stdev > (2.5 * float(dev))) and good:
                    print("\nFragment {}:".format(fragment))
                    print(
                            '*** Suspicious restraints in SADI line {} with high standard deviation {:4.3f} '
                            '(median length: {:4.3f} A) ***'.format(num + 1, stdev, median(distances)))
                    print('*** ' + restraints[num] + ' ***\n')
                    good = False
        if good:
            return True
        else:
            return False

    def search_for_error_response(self, fragment):
        """
        searches for a fragment name in the db as response to an invalid fragment name.
        :param fragment: the fragment
        :type fragment: string
        """
        result = search_fragment_name(fragment.lower(), self.databases)
        print('Do you mean one of these?:\n')
        print_search_results(result)

    def get_atoms(self, fragment, invert=False, cartesian=False):
        # type: (str, bool, bool) -> list
        """
        returns the atoms from the dbentry:
        [['O1', '1', '0.01453', '-1.6659', '-0.10966'],
        ['C1', '1', '0.00146', '-0.26814', '-0.06351'], ... ]
        :param fragment: fragment name
        :type fragment: string

        >>> db = ParseDB('./dsr_db.txt')
        >>> db.get_atoms('water', cartesian=False) # doctest: +NORMALIZE_WHITESPACE
        [['O1', 4, 0.0, 0.0, 0.0], ['H1', 2, 0.9584, 0.0, 0.0], ['H2', 2, -0.2392, 0.9281, 0.0]]
        >>> db.get_atoms('isoprop', cartesian=True) # doctest: +NORMALIZE_WHITESPACE
        [['O1', 1, 5.611626689944858, 4.650686219462244, 9.549527472636036],
        ['C1', 1, 5.969103603542006, 3.2932494073414893, 9.793246423730249],
        ['C2', 1, 4.814135175316713, 2.3627923165385942, 9.795180701119886],
        ['C3', 1, 7.0650724256725805, 2.9180586416449548, 8.82997628369121]]
        """
        fragment = fragment.lower()
        try:
            atoms = deepcopy(self.databases[fragment]['atoms'])
        except KeyError:
            print('*** Could not find {} in database ***'.format(fragment))
            sys.exit()
        if cartesian:
            cell = self.get_cell(fragment)
            for num, c in enumerate(self.get_coordinates(fragment)):
                atoms[num][2:5] = frac_to_cart(c, cell)
        if invert:
            atoms = invert_atomic_coordinates(atoms)
        return atoms

    def get_cell(self, fragment):
        # type: (str) -> list
        """
        returns the line with FRAG 17 cell from the dbentry
        >>> db = ParseDB('./dsr_db.txt')
        >>> db.get_cell('water')
        [1.0, 1.0, 1.0, 90.0, 90.0, 90.0]
        >>> db.get_cell('foobar')
        Traceback (most recent call last):
        ...
        SystemExit
        """
        try:
            cell = [float(x) for x in self.databases[fragment.lower()]['cell']]
        except KeyError:
            print(not_existing_error.format(fragment))
            sys.exit()
        return cell

    def get_startline(self, fragment):
        # type: (str) -> int
        """
        returns the line number from the dbentry
        """
        try:
            return int(self.databases[fragment.lower()]['startline'])
        except KeyError:
            print(not_existing_error.format(fragment))
            return 0
        except ValueError:
            print('Could not find startline of Fragment {}'.format(fragment))
            return 0

    def get_coordinates(self, fragment, cartesian=False, invert=False):
        # type: (str, bool, bool) -> list
        """
        Returns the coordinates of the fragment. Optionally is direct conversion to cartesian coordinates.
        [['0.01453', '-1.6659', '-0.10966'], [...] ]
        :param fragment:
        :param cartesian:
        :return:
        """
        try:
            atoms = self.get_atoms(fragment, invert=invert, cartesian=cartesian)
            coords = [i[2:5] for i in atoms]
        except KeyError:
            print(not_existing_error.format(fragment))
            sys.exit()
        return coords

    def get_restraints(self, fragment: str) -> list[str]:
        """
        returns the header of the dbentry of fragment.
        This header does not include comments, only the restraints.
        """
        fragment = fragment.lower()
        try:
            restr = [x.upper() for x in self.databases[fragment]['restraints']]
        except KeyError:
            print(not_existing_error.format(fragment))
            sys.exit()
        return restr

    def get_resi(self, fragment):
        """
        returns the residue name of the dbentry of fragment
        can be either class or class + number.
        convention is only class.
        """
        try:
            return self.databases[fragment.lower()]['resi']
        except KeyError:
            print(not_existing_error.format(fragment))
            return ''

    def get_fragment_name(self, fragment):
        """
        returns the first comment line of the dbentry of a fragment
        if a line with "rem Name:" is present, this line is used as comment.
        :param fragment: actual fragment name
        :type fragment: string
        """
        try:
            name = self.databases[fragment.lower()]['name']
        except KeyError:
            print(not_existing_error.format(fragment))
            return ''
        return name

    def get_src(self, fragment):
        """
        returns the source line of the dbentry of a fragment
        if a line with "rem Src:" is present.
        :param fragment: actual fragment name
        :type fragment: string
        """
        try:
            src = self.databases[fragment.lower()]['source']
        except KeyError:
            print(not_existing_error.format(fragment))
            return ''
        return src

    def get_db_name(self, fragment):
        """
        returns the fragment database name of fragment x
        """
        return self.databases[fragment.lower()]['dbname']

    def get_fragment_tags(self):
        tags = [x[0] for x in self.list_fragments()]
        return tags


def get_first_last_atom(atoms):
    """
    :type atoms: list
    returns the first and the last atom from the imported atom list
    """
    try:
        first = atoms[0][0]
        last = atoms[-1][0]
    except (KeyError, IndexError):
        return False
    return first, last


class ImportGRADE():
    def __init__(self, grade_tar_file, db, invert=False, maindb=None, userdb=None):
        """
        class to import fragments from GRADE of Global Phasing Ltd.
        :param grade_tar_file:
        :type grade_tar_file:
        :param invert:
        :type invert:
        :type maindb:   string
        :param userdb:  directory where the user database is located.
                        Default is the users home directory.
        :type userdb:   string
        """
        self.user_db_path = userdb
        self.main_db_path = maindb
        self.el = Element()
        self.invert = invert
        self._db = db
        self._db_dir = os.path.expanduser("~")
        self._db_tags = self._db.get_fragment_tags()
        self._db = self._db.databases
        gradefiles = self.get_gradefiles(grade_tar_file)
        self._pdbfile = gradefiles[0]
        self._dfixfile = gradefiles[1]
        # self._obpropfile = gradefiles[2]
        self._atoms = self.get_pdbatoms(self._pdbfile)
        self._firstlast = get_first_last_atom(self._atoms)
        self._restraints = self.get_restraints()
        self._resi_name = self.get_resi_from_pdbfile()
        if not isinstance(self._resi_name, str):
            self._resi_name = self._resi_name.decode()

    def get_gradefiles(self, grade_tar_file):
        """
        returns the .mol2 and .dfix file location
        """
        grade_base_filename = os.path.splitext(grade_tar_file)
        if grade_base_filename[1] == '.tgz':
            try:
                gradefile = tarfile.open(grade_tar_file)  # , encoding="ascii")
            except IOError:
                print('No such file or directory: {}'.format(grade_tar_file))
                sys.exit()
        else:
            print('*** File {} is not a valid file to import from ***'.format(grade_base_filename[1]))
            sys.exit()
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

    @staticmethod
    def _read_tarfile_as_ascii(filehandle):
        tmp = []
        for line in filehandle:
            if not isinstance(line, str):
                line = line.decode('ascii')
            tmp.append(line)
        return tmp

    def get_name_from_pdbfile(self):
        # type: () -> str
        """
        get the fragment name from the pdbfile.txt file

        >>> db = ParseDB('./dsr_db.txt')
        >>> mog = ImportGRADE('./tests/test-data/ALA.gradeserver_all.tgz', db)
        >>> mog.get_name_from_pdbfile()
        'Alanine'
        """
        full_name = None
        full_name_regex = re.compile(r'^.*Compound full name.*')
        for line in self._pdbfile:
            if not isinstance(line, str):
                line = line.decode('ascii')
            if not line:
                continue
            if full_name_regex.match(line):
                line = line.replace('_', '')
                line = line.replace('-', '')
                line = line.replace('#', '')
                try:
                    line = line.split()
                    _ = line[4]
                except IndexError:
                    full_name = ''
                else:
                    full_name = line[5]
        return full_name

    def get_resi_from_pdbfile(self):
        """
        get the fragment name from the pdbfile.txt file

        >>> db = ParseDB('./dsr_db.txt')
        >>> mog = ImportGRADE('./tests/test-data/ALA.gradeserver_all.tgz', db)
        >>> mog.get_resi_from_pdbfile()
        'ALA'
        """
        resi_name = None
        resi_regex = re.compile(r'^HETATM\s+1.*')
        for line in self._pdbfile:
            if not isinstance(line, str):
                line = line.decode('ascii')
            if not line:
                continue
            if resi_regex.match(line):
                line = line.replace('_', '')
                line = line.replace('-', '')
                line = line.replace('#', '')
                resi_name = line.split()[3]
                break
        return resi_name

    def get_comments(self):
        """
        returns a detailed comment out of the .dfix file where the fragment was imported from

        #   REM Produced by Grade Web Server http://grade.globalphasing.org
        #   REM GEN: Generated by GRADE 1.2.5 (December 20 2013)
        #   REM GEN: from SMILES C1=CC=C2C(=C1)C=C(C3=CC=CC=C23)Br
        #   REM GEN: using MOGUL 1.6(RC5), CSD as535be with quantum mechanics RM1
        #   REM grade-cif2shelx output
        #   REM Version: 0.0.5 <Dec 20 2013>
        """
        matches = ['REM Produced by Grade', 'REM GEN:', 'REM grade-cif2shelx', 'REM Version:',
                   'REM Total charge', 'REM Name:']
        comments = []
        name = 'REM Name: ' + self.get_name_from_pdbfile()
        comments.append(name.split())
        for m in matches:
            for line in self._dfixfile:
                if not isinstance(line, str):
                    line = line.decode()
                if line.startswith(m):
                    comments.append(line.split())
        return comments

    def get_restraints(self):
        """
        reads the restraints from a dfix-file
        """
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
        """
        returns the atoms from a pdb file
        """
        atomlines = []
        for line in pdb:
            if not isinstance(line, str):
                line = line.decode()
            line = line.split()
            if line[0] == 'HETATM':
                if line[-1] == 'H':
                    continue
                tmp = [line[2], line[11]]
                tmp.extend(line[6:9])
                atomlines.append(tmp)
        if not atomlines:
            self.import_error(pdb)
        if self.invert:
            atomlines = invert_atomic_coordinates(atomlines)
        return atomlines

    def bild_grade_db_entry(self):
        """
        builds a dbentry from the information supplied by GRADEs
        .mol2 and .dfix file
        """
        db_import_dict = {}
        num = 0
        name = self._resi_name[:4].upper()
        if not isinstance(name, str):
            name = name.decode()
        resi_name = name
        if not self._db_tags:
            print('*** Unable to import fragment. Database is empty. ***')
            sys.exit()
        # Check if anything is already in database:
        # tags are always lower case:
        while resi_name.lower() in self._db_tags:
            # Already there, so add a number to the tag:
            num = num + 1
            resi_name = resi_name[:3] + str(num)
        if num == 0:
            resi_name = resi_name[:3]
        else:
            resi_name = resi_name[:3] + str(num)
        fragline = 'FRAG 17 1  1  1  90  90  90'
        db_import_dict[resi_name] = {
            'restraints': self._restraints,
            'resi'      : resi_name,
            'cell'      : fragline.split(),
            'atoms'     : self._atoms,
            'line'      : None,
            'db'        : 'dsr_user_db',
            'comments'  : self.get_comments(),
            'name'      : resi_name
        }
        return db_import_dict

    @staticmethod
    def import_error(filename):
        # type: (str) -> NotImplemented
        """
        warns for import errors
        """
        print('*** Unable to import GRADE file {} GRADE import to DSR '
              'relies on GRADE v1.100 and up. ***'.format(filename))
        sys.exit()

    def write_user_database(self):
        """
        writes content of existing dsr_user_db.txt and the imported GRADE entry to
        the dsr_user_db.txt
        """
        fragments = list(self._db)
        imported_entry = self.bild_grade_db_entry()
        grade_db_names = list(imported_entry)
        user_db_names = []
        for i in fragments:
            if self._db[i]['dbname'] == 'dsr_user_db':
                user_db_names.append(self._db[i]['name'])
        # try to write the new dbentry:
        try:
            with open(self.user_db_path, 'w') as f:
                for name in grade_db_names:
                    print('Importing {} ({}) to user database...'.format(self._resi_name, name))
                    atomlist = imported_entry[name]['atoms']
                    comments = imported_entry[name]['comments']
                    comments = '\n'.join([' '.join(i) for i in comments if i])
                    head = '\n'.join([' '.join(x) for x in imported_entry[name]['restraints']])
                    atoms = '\n'.join(['{:<6} -{}  {:>8.3f}{:>8.3f}{:>8.3f}' \
                                      .format(y[0], self.el.get_atomic_number(y[1]), float(y[2]), float(y[3]),
                                              float(y[4])) for y in atomlist])
                    resi_name = str(name)
                    cell = '  '.join(imported_entry[name]['cell'])
                    dbentry = '<{}> \n{} \nRESI {} \n{} \n{} \n{} \n</{}>\n''\
                        '.format(resi_name, comments, resi_name, head, cell, atoms, resi_name)
                    f.write(dbentry)
        except IOError as e:
            print(e)
            sys.exit()
        # try to write existing dbentries:
        try:
            with open(self.user_db_path, 'a+') as fu:
                for tag in fragments:
                    if self._db[tag]['dbname'] == 'dsr_user_db':
                        atomlist = self._db[tag]['atoms']
                        head = '\n'.join([''.join(x) for x in self._db[tag]['restraints']])
                        atoms = '\n'.join(['{:<6}{:<2}{:>8.3f}{:>8.3f}{:>8.3f}'
                                          .format(y[0], y[1], float(y[2]), float(y[3]), float(y[4])) for y in atomlist])
                        resi_name = self._db[tag]['resi']
                        comments = '\n'.join(self._db[tag]['comments'])
                        name = self._db[tag]['name']
                        comments = "REM Name: {}\n{}".format(name, comments)
                        fragline = 'FRAG 17 {} {} {} {} {} {}'.format(*self._db[tag]['cell'])
                        dbentry = '\n<{}> \n{} \nRESI {} \n{} \n{} \n{} \n</{}>\n' \
                                  ''.format(tag, comments, resi_name, head, fragline, atoms, tag)
                        fu.write(dbentry)
        except IOError as e:
            print(e)
            sys.exit()
        print('User database successfully updated.')


if __name__ == '__main__':
    """
    import doctest
    failed, attempted = doctest.testmod()  # verbose=True)
    if failed == 0:
        print('passed all {} tests!'.format(attempted))
    else:
        print('{} of {} tests failed'.format(failed, attempted))
    #homedir = expanduser("~")
    #userdb_path = os.path.join(homedir, "dsr_db.txt")
    """


    def std_dev_test(fragment, restraints, atoms, atnames, cell):
        for num, line in enumerate(restraints):
            prefixes = []
            dev = 0.02
            line = line.split()
            if not line:
                continue
            if line[0].upper() == 'SADI':
                prefixes.append(line[0])
                del line[0]
                try:
                    if not str(line[0][0]).isalpha():
                        prefixes.append(line[0])
                        dev = line[0]
                        del line[0]  # delete standard deviation
                except IndexError:
                    return False
                if len(line) % 2 == 1:  # test for uneven atoms count
                    print('*** Inconsistent SADI restraint line {} of "{}". '
                          'Not all atoms form a pair ***'.format(num, fragment))
                pairs = pairwise(line)
                distances = []
                pairlist = []
                if len(pairs) <= 2:
                    return True
                for i in pairs:
                    pairlist.append(i)
                    try:
                        atom1 = [float(x) for x in atoms[atnames.index(i[0])][2:5]]
                        atom2 = [float(x) for x in atoms[atnames.index(i[1])][2:5]]
                    except ValueError:
                        return False
                    dist = atomic_distance(atom1, atom2, cell)
                    distances.append(dist)
                if len(pairlist) == 2:
                    pairdev = 1 - (min(distances) / max(distances))
                    if pairdev > 2 * dev:
                        print('foooooooo', pairdev)
                stdev = std_dev(distances)
                print("esd < 0.065 ?: {:<4.4f}, esd < 2.5*sigma?: {:<4.4f} < {:<4.4f} -> {}".format(
                        stdev, stdev, 2.5 * float(dev), ("yes" if stdev < 2.5 * float(dev) else "no")))


    frag = 'WBVNT'

    dbpath = os.path.abspath('./dsr_db.txt')
    userpath = os.path.abspath('c:/Users/daniel/dsr_user_db.txt')
    db = ParseDB(dbpath, userdb_path=userpath)
    # pprint(db.databases['toluene'])
    db.check_consistency(frag)
    db.check_db_atom_consistency(frag)
    db.check_db_restraints_consistency(frag)
    db.check_sadi_consistence(frag)
    atnames = db.get_atomnames(frag)
    restr = db.get_restraints(frag)
    atoms = db.get_atoms(frag)
    cell = db.get_cell(frag)
    # print(db.databases['toluene'])
    print(dbpath)

    std_dev_test(fragment=frag, restraints=restr, atnames=atnames, atoms=atoms, cell=cell)
