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

import os
import re
import sys
import tarfile
from collections import Counter
from copy import deepcopy
from os.path import expanduser
from pprint import pprint

import misc
from atomhandling import get_atomtypes
from atoms import Element
from constants import atomregex, SHX_CARDS, RESTRAINT_CARDS, sep_line
from misc import atomic_distance, nalimov_test, std_dev, median, pairwise, \
    unwrap_head_lines, dice_coefficient2
from misc import touch

atreg = re.compile(atomregex)


def invert_dbatoms_coordinates(atoms):
    """
    Inverts SHELXL atom coordinates like
    [[C1  1  0.44  0.21  -1.23 ][ ...]]

    :param atoms: list of atom list
    """
    for line in atoms:
        try:
            inv_coord = [str(-float(i)) for i in line[-3:]]
        except:
            print('Unable to invert fragment coordinates.')
            return False
        line[-3:] = inv_coord
    return atoms


def search_fragment_name(search_string, db, numresults=6):
    """
    searches the Name: comments in the database for a given name

    :param numresults: number of results to return after search
    :type search_string: str
    :param gdb: database object
    :type gdb: global_DB
    """
    names_list = []
    for fragment in db:
        fragname = db[fragment]['name']
        line_number = db[fragment]['startline']
        dbname = db[fragment]['dbname']
        names_list.append([fragment, fragname, line_number, dbname])
    search_results = []
    for fragment in names_list:
        key = make_sortkey(fragment[1], searchkey=True)
        # key[0] and key[1] combined to have also the sum formula:
        coefficient = dice_coefficient2(search_string, "{0}{1}".format(key[0], key[1]))
        fragment.append([coefficient, key[1]])
        search_results.append(fragment)
    # select the best n results, sort for first and sewcond search key
    selected_results = sorted(search_results, key=lambda coeff: [coeff[-1][0], coeff[-1][1]], reverse=False)[
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


__metaclass__ = type  # use new-style classes


# hardwired names of the database files:
# dsr_db.txt is the db from the distribution. This file should not be edited.
# dsr_user_db.txt if this file exists, all its content is also read in.

def read_file_data(filepath):
    # type: (str, bool) -> list
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
        # type: () -> NotImplemented
        """
        self._databases: dictionary with the individial fragments
        """
        self.databases = {}
        self.maindb_path = maindb_path
        self.userdb_path = userdb_path
        if maindb_path:
            self.parse(maindb_path, 'dsr_db')
        if userdb_path:
            self.parse(userdb_path, 'dsr_usr_db')

    def parse(self, dbpath, dbname):
        # type: (str, str) -> dict
        """
        This method returns all fragment name tags in the database

        >>> dbpath = os.path.abspath('dsr_db.txt')
        >>> db = ParseDB(dbpath)
        >>> db.parse('dsr_db')['water']
        ... # doctest: +NORMALIZE_WHITESPACE
        {'startline': 2402,
        'name': 'Water, H2O',
        'comments': [],
        'atoms': [['O1', '4', '0.0000', '0.0000', '0.0000'],
                  ['H1', '2', '0.9584', '0.0000', '0.0000'],
                  ['H2', '2', '-0.2392', '0.9281', '0.0000']],
        'cell': [1.0, 1.0, 1.0, 90.0, 90.0, 90.0],
        'source': 'pbe1pbe/6-311++G(3df,3pd), Ilia A. Guzei',
        'resi': 'H2O',
        'restraints': ['DFIX 0.9584 0.001 O1 H1 O1 H2',
                       'DFIX 1.5150 0.001 H1 H2'],
        'endline': 2412,
        'dbname': 'dsr_db'}
        """
        frag = ''
        db = {}
        start_regex = re.compile(r'<[^/].*>', re.IGNORECASE)  # regular expression for db tag.
        starttag = False
        fraglines =[]
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
                db[frag].update(
                    {'endline': num,
                     'startline': startnum})
                db = self.parse_fraglines(frag, fraglines, db)
                fraglines = []
            # start tag was found, appending lines to fragment list
            if starttag:
                fraglines.append(line)
            # matching start tag and compiling end regex
            if start_regex.match(line):
                if starttag:
                    print('*** Error in database "{}.txt" in line {}. End tag is missing ***'.format(dbname, num+1))
                frag = line.strip('<> \n\r').lower()
                if frag in db:
                    print('\n*** Duplicate database entry "{}" found! Please remove/rename '
                          'second entry\nand/or check all end tags in the database dsr_usr_db.txt '
                          'or dsr_db.txt. ***'.format(frag))
                starttag = True
                startnum = num
                db[frag] = {'dbname': dbname}
                end_regex = re.compile(re.escape(r'</{}>'.format(frag)), re.IGNORECASE)
                continue
        self.databases.update(db)
        return db

    @staticmethod
    def parse_fraglines(fragname, fraglines, db):
        # type: (str, list, dict) -> dict
        """
        Fills the database dictionary with fragment data.
        """
        fragname = fragname.lower()
        headlist = []
        comments = []
        residue = ''
        atoms = []
        cell = []
        source = ''
        name = ''
        nameregex = re.compile(r'REM\s+NAME:', re.IGNORECASE)
        rem = re.compile(r"REM.*", re.IGNORECASE)
        srcregex = re.compile(r'REM\s+(SRC:|SOURCE:)', re.IGNORECASE)
        # devide atoms and the rest:
        for num, aline in enumerate(fraglines):
            if atreg.match(aline):  # search atoms
                atline = aline.split()[:5]  # convert to list and use only first 5 columns
                if atline[0] not in SHX_CARDS:  # exclude all non-atom cards
                    atoms.append(atline)
                    fraglines[num] = ''
        # bring all wrapped lines with = at end to a single line:
        fraglines = unwrap_head_lines(fraglines)
        if fragname == 'cpstar':
            pass
        for num, line in enumerate(fraglines):
            # collect the comments:
            if rem.match(line):
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
            com = line[:4].upper()
            if com and com in SHX_CARDS:
                if line[:4].upper() == 'RESI':  # faster than startswith()?
                    resiline = line.split()
                    for n, i in enumerate(resiline):
                        if not i[0].isalpha():
                            del resiline[n]
                    try:
                        residue = resiline[1].upper()
                        continue
                    except IndexError:
                        print('*** Invalid residue definition in database entry {}. ***'.format(fragname))
                        sys.exit()
                # collect the unit cell of the fragment:
                if line[:4].upper().startswith('FRAG'):
                    try:
                        cell = [float(x) for x in line.split()[2:8]]
                    except ValueError:
                        print('*** Invalid unit cell found in database entry {}. ***'.format(fragname))
                    continue
                # these must be restraints:
                if line[:4] in RESTRAINT_CARDS:
                    headlist.append(line)
            elif com:
                print('*** Bad line {} in header of database entry "{}" found! ({}.txt) ***'
                      .format(num+db[fragname]['startline']+2, fragname, db[fragname]['dbname']))
                print(line)
        db[fragname].update({
            'restraints': headlist,  # header with just the restraints
            'resi': residue,  # the residue class
            'cell': cell,  # FRAG ...
            'atoms': atoms,  # the atoms as lists of list
            'comments': comments,  # the comment line
            'source': source,
            'name': name})
        if not db:
            print('*** No database found! ***\n')
            sys.exit()
        return db

    def __getitem__(self, fragment):
        try:
            return self.databases[fragment]
        except KeyError:
            print("*** Fragment {} was not found in the Database! ***".format(fragment))
            sys.exit()

    def list_fragments(self):
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
        """
        returns the atomic numbers of the atoms in a fragment in
        same order as the dsr db as list
        :param fragment: fragment name tag
        :type fragment: basestring
        """
        atnumbers = []
        types = get_atomtypes(self.get_atoms(fragment))
        for i in types:
            el = Element()
            atnumbers.append(el.get_atomic_number(i))
        return atnumbers

    def get_atomnames(self, fragment=''):
        atoms = self.databases[fragment]['atoms']
        names = []
        for i in atoms:
            names.append(i[0])
        return names

    def get_head_for_gui(self, fragment):
        """
        returns header information of the specific fragment:
        tag, Name/comment, source, cell, residue, dbtype, restr, atoms
        >>> dbpath = os.path.abspath('dsr_db.txt')
        >>> db = ParseDB()
        >>> p = db.parse(dbpath, 'dsr_db')
        >>> db.get_head_for_gui('benZene')
        ... # doctest: +NORMALIZE_WHITESPACE
        <tag>
         benzene
        </tag>
        <comment>
         Benzene, Benzol, C6H6
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
         SADI 0.02 C1 C2 C2 C3 C3 C4 C4 C5 C5 C6 C6 C1;;SADI 0.04 C1 C5 C1 C5 C4 C2 C4 C6 C6 C2 C5 C3;;FLAT C1 > C6;;SIMU C1 > C6;;RIGU C1 > C6
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
        self.check_db_header_consistency(fragment)
        if not self.check_sadi_consistence(fragment):
            sys.exit()
        self.check_db_atom_consistency(fragment)

    def check_consistency(self, fragment):
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
        if not dbentry['restraints']:
            print('*** Restraints in the header of database entry "{}" ({}) missing! Check your Database! ***' \
                  .format(fragment, dbentry['name']))
            return False
        if not dbentry['resi']:
            print('*** Residue in the header of database entry "{}" ({}) missing! Check your Database! ***' \
                  .format(fragment, dbentry['name']))
        if not dbentry['atoms']:
            print('*** Database entry of "{}" in line {} of "{}.txt" is corrupt. '
                  'No atoms found! ***'.format(fragment, dbentry['startline']+1, dbentry['dbname']))
            print('*** Have you really followed the syntax? ***')
            sys.exit()
        if not dbentry['endline']:
            print('*** Could not find end of dbentry for fragment "{}" in line {} of "{}.txt". '
                  'Check your database. ***'.format(fragment, dbentry['startline']+1, dbentry['dbname']))
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
        """
        This method is for atoms only (without db header)!
        check the db for duplicates:
        """
        fragment = fragment.lower()
        dbatoms = [i[0].upper() for i in self.databases[fragment]['atoms']]
        # check for duplicates:
        while dbatoms:
            at = dbatoms.pop()
            if at in dbatoms:
                print('*** Duplicate atom {0} in database entry "{1}" ({2}) '
                      'found! Check your database... ***'
                      .format(at, fragment, self.databases[fragment]['name']))
                sys.exit()

    def check_db_header_consistency(self, fragment):
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
                if i in SHX_CARDS:
                    continue
                if i in ('>', '<'):
                    continue
                try:
                    float(i)
                except ValueError:
                    restraint_atoms_list.add(i)
        for atom in restraint_atoms_list:
            if atom.upper() not in atoms:
                status = False
                print('\n*** Unknown atom "{}" in restraints of "{}: {}" in database {}. ***'
                      .format(atom, fragment, self.get_fragment_name(fragment), self.get_db_name(fragment)))
        if not status:
            print('*** Check database entry. ***\n')
            #sys.exit()
        return status

    def check_sadi_consistence(self, fragment):
        """
        check if same distance restraints make sense. Each length of an atom
        pair is tested agains the standard deviation of all distances.
        For a large standard deviation, the list is tested for outliers.
        :param atoms: atoms list of thr fragment
        :param restraints: restraints list
        :param fragment: frag name
        :param factor: factor for confidence interval
        """
        fragment = fragment.lower()
        atoms = self.get_atoms(fragment)
        restr = self.get_restraints(fragment)
        restraints = deepcopy(restr)
        atnames = [i[0].upper() for i in atoms]
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
                        a = atoms[atnames.index(i[0])][2:5]
                        b = atoms[atnames.index(i[1])][2:5]
                    except ValueError:
                        return False
                    a = [float(x) for x in a]
                    b = [float(y) for y in b]
                    dist = atomic_distance(a, b, self.get_cell(fragment))
                    distances.append(dist)
                stdev = std_dev(distances)
                # only do outlier test if standard deviation is suspiciously large:
                if stdev > 0.065:
                    outliers = nalimov_test(distances)
                    if outliers:
                        print("\nFragment {}:".format(fragment))
                        for x in outliers:
                            pair = ' '.join(pairlist[x])
                            print(
                                '*** Suspicious deviation in atom pair "{}" ({:4.3f} A, median: {:4.3f}) of '
                                'SADI line {} ***'.format(pair, distances[x], median(distances), num + 1))
                            print('*** ' + restr[num][:60] + ' ... ***')
                            return False
                if stdev > 2.5 * float(dev):
                    print("\nFragment {}:".format(fragment))
                    print(
                        '*** Suspicious restraints in SADI line {} with high standard deviation {:4.3f} '
                        '(median length: {:4.3f} A) ***'.format(num + 1, stdev, median(distances)))
                    print('*** ' + ' '.join(prefixes + line) + ' ***')
                    return False
        return True

    def search_for_error_response(self, fragment):
        """
        searches for a fragment name in the db as response to an invalid fragment name.
        :param fragment: the fragment
        :type fragment: string
        """
        result = search_fragment_name(fragment.lower(), self.databases)
        print('Do you mean one of these?:\n')
        print_search_results(result)

    def get_atoms(self, fragment):
        """
        returns the atoms from the dbentry:
        [['O1', '1', '0.01453', '-1.6659', '-0.10966'],
        ['C1', '1', '0.00146', '-0.26814', '-0.06351'], ... ]
        :param fragment: fragment name
        :type fragment: string
        """
        fragment = fragment.lower()
        try:
            return self.databases[fragment]['atoms']
        except KeyError:
            print('*** Could not find {} in database ***'.format(fragment))
            sys.exit()

    def get_cell(self, fragment):
        """
        returns the line with FRAG 17 cell from the dbentry
        """
        fragment = fragment.lower()
        try:
            cell = self.databases[fragment]['cell']
        except KeyError:
            print('*** Fragment "{}" not found in database ***'.format(fragment))
            sys.exit()
        return cell

    def get_startline(self, fragment):
        """
        returns the line number from the dbentry
        """
        return self.databases[fragment.lower()]['startline']

    def get_restraints(self, fragment):
        """
        returns the header of the dbentry of fragment.
        This header does not include comments, only the restraints.

        Parameters
        ----------
        fragment : str

        Returns
        ----------
        list
        """
        head = []
        fragment = fragment.lower()
        try:
            head = self.databases[fragment]['restraints']
        except KeyError:
            print('*** Could not find {} in database ***'.format(fragment))
            sys.exit()
        return head

    def get_resi(self, fragment):
        """
        returns the residue name of the dbentry of fragment
        can be either class or class + number.
        convention is only class.
        """
        return self.databases[fragment.lower()]['resi']

    def get_fragment_name(self, fragment):
        """
        returns the first comment line of the dbentry of a fragment
        if a line with "rem Name:" is present, this line is used as comment.
        :param fragment: actual fragment name
        :type fragment: string
        """
        name = self.databases[fragment.lower()]['name']
        return name

    def get_src(self, fragment):
        """
        returns the source line of the dbentry of a fragment
        if a line with "rem Src:" is present.
        :param fragment: actual fragment name
        :type fragment: string
        """
        src = self.databases[fragment.lower()]['source']
        return src

    def get_db_name(self, fragment):
        """
        returns the fragment database name of fragment x
        """
        return self.databases[fragment.lower()]['dbname']


class ImportGRADE():
    def __init__(self, grade_tar_file, invert=False, maindb=None, userdb=None):
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
        :param dbnames: file names of the databases
        :type dbnames:  string
        """
        if userdb is None:
            homedir = expanduser("~")
            self.user_db_path = os.path.join(homedir, "dsr_user_db.txt")
            if not os.path.isfile(self.user_db_path):
                touch(self.user_db_path)
        else:
            self.user_db_path = userdb
        ##############################################    
        if maindb is None:
            try:
                main_dbdir = os.environ["DSR_DIR"]
            except(KeyError):
                main_dbdir = './'
            self.main_db_path = os.path.join(main_dbdir, 'dsr_db.txt')
        else:
            self.main_db_path = maindb
        ##############################################        
        self.el = Element()
        self.invert = invert
        self._getdb = ReadDB(self.main_db_path, self.user_db_path)
        self._db_dir = expanduser("~")
        self._db_tags = self._getdb.find_db_tags()
        self._gdb = global_DB(invert=False)
        self._db = self._gdb.build_db_dict()
        gradefiles = self.get_gradefiles(grade_tar_file)
        self._pdbfile = gradefiles[0]
        self._dfixfile = gradefiles[1]
        # self._obpropfile = gradefiles[2]
        self._atoms = self.get_pdbatoms(self._pdbfile)
        self._firstlast = self.get_first_last_atom(self._atoms)
        self._restraints = self.get_restraints()
        self._resi_name = self.get_resi_from_pdbfile()
        self._frag_name = self.get_name_from_pdbfile()
        if not isinstance(self._resi_name, str):
            self._resi_name = self._resi_name.decode()
        self._comments = self.get_comments()

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
        """
        get the fragment name from the pdbfile.txt file
        :param pdbfile: file with some information about the molecule
        :type pdbfile: list of strings

        >>> mog = ImportGRADE('./test-data/ALA.gradeserver_all.tgz')
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
                except(IndexError):
                    full_name = None
                full_name = line[5]
                break
        return full_name

    def get_resi_from_pdbfile(self):
        """
        get the fragment name from the pdbfile.txt file
        :param pdbfile: file with some information about the molecule
        :type pdbfile: list of strings

        >>> mog = ImportGRADE('./test-data/ALA.gradeserver_all.tgz')
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
        matches = ['REM Produced by Grade', 'REM GEN:', 'REM grade-cif2shelx', 'REM Version:', 'REM Total charge']
        comments = []
        name = 'REM Name: ' + self._frag_name
        comments.append(name.split())
        for m in matches:
            for line in self._dfixfile:
                if not isinstance(line, str):
                    line = line.decode()
                if line.startswith(m):
                    comments.append(line.split())
        return comments

    def get_first_last_atom(self, atoms):
        """
        :type atoms: list
        returns the first and the last atom from the imported atom list
        """
        try:
            first = atoms[0][0]
            last = atoms[-1][0]
        except:
            return False
        return first, last

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
            atomlines = invert_dbatoms_coordinates(atomlines)
        return atomlines

    def bild_grade_db_entry(self):
        """
        builds a dbentry from the information supplied by GRADEs
        .mol2 and .dfix file
        """
        db_import_dict = {}
        num = 0
        name = self._resi_name[:3].upper()
        if not isinstance(name, str):
            name = name.decode()
        resi_name = name
        if not self._db_tags:
            print('*** Unable to import fragment. Database is empty. ***')
            sys.exit()
        for i in self._db_tags:
            while resi_name == i[0]:
                num = num + 1
                resi_name = resi_name[:3] + str(num)
        if num == 0:
            resi_name = resi_name[:3]
        else:
            resi_name = resi_name[:3] + str(num)
        fragline = 'FRAG 17 1  1  1  90  90  90'
        db_import_dict[resi_name] = {
            'head': self._restraints,
            'resi': resi_name,
            'fragline': fragline.split(),
            'atoms': self._atoms,
            'line': None,
            'db': 'dsr-user-db',
            'comment': self._comments,
            'name': resi_name
        }
        return db_import_dict

    def import_error(self, filename):
        """
        warns for import errors
        """
        print('*** Unable to import GRADE file {}'.format(filename))
        print('GRADE import relies on GRADE v1.100 and up. ***')
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
            if self._db[i]['db'] == 'dsr_user_db':
                user_db_names.append(self._db[i]['name'])
        # try to write the new dbentry:
        try:
            with open(self.user_db_path, 'w') as f:
                for name in grade_db_names:
                    print('Importing {} ({}) to user database...'.format(self._resi_name, name))
                    atomlist = imported_entry[name]['atoms']
                    comment = imported_entry[name]['comment']
                    comment = '\n'.join([' '.join(i) for i in comment if i])
                    head = '\n'.join([' '.join(x) for x in imported_entry[name]['head']])
                    atoms = '\n'.join(['{:<6} -{}  {:>8.3f}{:>8.3f}{:>8.3f}' \
                                      .format(y[0], self.el.get_atomic_number(y[1]), float(y[2]), float(y[3]),
                                              float(y[4])) for y in atomlist])
                    resi_name = str(name)
                    cell = '  '.join(imported_entry[name]['fragline'])
                    dbentry = '<{}> \n{} \nRESI {} \n{} \n{} \n{} \n</{}>\n''\
                        '.format(resi_name, comment, resi_name, head, cell, atoms, resi_name)
                    f.write(dbentry)
        except IOError as e:
            print(e)
            sys.exit()
        # try to write existing dbentries:
        try:
            with open(self.user_db_path, 'a+') as fu:
                for i in fragments:
                    name = i
                    if self._db[i]['db'] == 'dsr_user_db':
                        # userdb = list(self._db[i])
                        atomlist = self._db[i]['atoms']
                        head = '\n'.join([''.join(x) for x in self._db[i]['head']])
                        atoms = '\n'.join(['{:<6}{:<2}{:>8.3f}{:>8.3f}{:>8.3f}' \
                                          .format(y[0], y[1], float(y[2]), float(y[3]), float(y[4])) for y in atomlist])
                        resi_name = self._db[i]['resi']
                        comment = self._db[i]['comment']
                        comment = '\nREM '.join(comment)
                        fragline = '  '.join(self._db[i]['fragline'])
                        dbentry = '\n<{}> \nREM {} \nRESI {} \n{} \n{} \n{} \n</{}>\n' \
                                  ''.format(name, comment, resi_name, head, fragline, atoms, name)
                        fu.write(dbentry)
        except IOError as e:
            print(e)
            sys.exit()
        print('User database successfully updated.')


if __name__ == '__main__':
    import doctest
    failed, attempted = doctest.testmod()  # verbose=True)
    if failed == 0:
        print('passed all {} tests!'.format(attempted))
    else:
        print('{} of {} tests failed'.format(failed, attempted))
    #homedir = expanduser("~")
    #userdb_path = os.path.join(homedir, "dsr_db.txt")
    dbpath = os.path.abspath('./dsr_db.txt')
    db = ParseDB(dbpath)
    pprint(db.databases['cpstar'])
    print(dbpath)