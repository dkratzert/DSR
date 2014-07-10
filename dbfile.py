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

from constants import atomregex, SHX_CARDS, RESTRAINT_CARDS
import misc



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


__metaclass__ = type  # use new-style classes

# hardwired names of the database files:
# dsr_db.txt is the db from the distribution. This file should not be edited.
# dsr_user_db.txt if this file exists, all its content is also read in.

class ReadDB():
    '''
    reads in the system and user db files (self._db_file_names) and makes
    a dictionary of them.
    '''

    def __init__(self, dbdir='', dbnames=''):
        self._db_file_names = dbnames
        try:
            self._db_dir = dbdir
        except(KeyError):
            print('\nThe environment variable DSR_DB_DIR was not found.\n'\
                'Please set this variable to the path of the DSR database!\n')
            self._db_dir = './'
        self._databases = self.getDB_files_dict()


    def getDBpath(self, db_file_name):
        '''
        returns the full db path as tuple
        '''
        fullpath = os.path.join(self._db_dir, db_file_name)  # full path with filename
        return fullpath


    def getDB_files_dict(self):
        '''
        returns the database as dictionary. Each file has its own key.
        {'dsr-db': ('line1\n', 'line2\n', '...'), 'dsr-user-db': ('line1\n', 'line2\n', '...')}
        '''
        db_dict = {}
        for name in self._db_file_names:
            dblist = []
            filename = self.getDBpath(name)
            base_filename = os.path.splitext(name)[0]
            try:
                with open(filename, 'r') as f:
                    for line in f:
                        dblist.append(line)
            except(IOError) as e:
                if not str(e).find('dsr_db'):
                    print(e)
                    sys.exit(-1)
                else:
                    continue
            db_dict[base_filename] = dblist
        return db_dict


    def find_db_tags(self):
        '''
        This method lists all fragment names in the database
        '''
        regex = r'^<[^/].*>'  # regular expression for db tag.
        dbnames = []
        dbkeys = list(self._databases.keys())  # names of the databases
        for db in dbkeys:
            for num, line in enumerate(self._databases[db]):
                if re.match(regex, line):
                    frag = [str(line.strip('<> \n\r')).upper(), num + 1, db]  # stripping spaces is important here!
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
                    'second entry\nand/or check all end tags.\n'.format(duplicates.pop()))
            sys.exit(False)
        # # sort lower-case:
        dbnames.sort(key=lambda x: x[0].lower())
        return dbnames



class global_DB():
    '''
    creates a final dictionary where all dbs are included:

    {'toluene':
       {'head'   : ('rem Source: Gaussian', 'SADI C1 C2 ...')
        'resi'   : 'TOL'
        'frag'   : 'FRAG 17 1 1 1 90 90 90'
        'atoms'  : ('C1 1 0.234234 1.567234 0.845662', 'C2 1 ...')
        'line'   : '123'
        'db'     : 'dsr_db'
        'comment': 'rem Source: Gaussian'
       }
     'benzene':
       {'head'   : ()
        'resi'   : ''
        'frag'   : ''
        'atoms   : ()
        'line'   : ''
        'db'     : 'dsr_db' or 'dsr_user_db'
        'comment': ''
       }
    }
    '''

    def __init__(self, invert=False, dbdir=os.environ["DSR_DB_DIR"],
                 dbnames=["dsr_db.txt", "dsr_user_db.txt"]):
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
        :param dbdir:  directory where the database is located. Default is environment variable DSR_DB_DIR
        :type dbdir:   string
        :param dbnames: file names of the databases
        :type dbnames:  string
        '''
        self.invert = invert
        self._getdb = ReadDB(dbdir, dbnames)
        self._db_tags = self._getdb.find_db_tags()
        self._db_plain_dict = self._getdb.getDB_files_dict()
        self._dbentry_dict = self.build_db_dict()


    def build_db_dict(self):
        '''
        returns the global db dictionary
        '''
        db_dict = {}
        for i in self._db_tags:
            fragment = i[0].lower()
            line = i[1]
            db = str(i[2])
            headboth = self.get_head_lines(fragment, db, line)
            header = headboth[0]
            fragline = headboth[1]
            residue = self.get_residue_from_head(header, fragment)
            atoms = self.get_fragment_atoms(fragment, db, line)
            comment = self.get_head_lines(fragment, db, line)[2]
            comment = [' '.join(x) for x in comment]
            db_dict[fragment] = {
                'head'    : header,
                'resi'    : residue,
                'fragline': fragline.split(),
                'atoms'   : atoms,
                'line'    : line,
                'db'      : db,
                'comment' : comment,
                'name'    : i[0]
                }
        return db_dict


    def get_residue_from_head(self, head, fragment=''):
        '''
        extracts the residue from the head
        and returns just the first string after RESI
        :param head: ['line1\n', 'line2\n', '...']
        :type head: list of strings
        '''
        for index, line in enumerate(head):
            if line.upper().startswith('RESI'):
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
        for i in self._db_plain_dict[db_name][int(line):]:
            i = i.strip('\n\r')
            regex = re.escape(r'</{}>'.format(fragment.lower()))
            if re.match(regex, i.lower()):  # find the endtag of db entry
                end = True
                break
            if re.search(atomregex, str(i)):  # search atoms
                l = i.split()[:5]  # convert to list and use only first 5 columns
                if l[0].upper() not in SHX_CARDS:  # exclude all non-atom cards
                    atoms.append(l)
        if not atoms:
            print('database entry of "{}" in line {} of "{}.txt" is corrupt. No atoms found!'.format(fragment, line, db_name))
            print('Have you really followed the syntax?')
            sys.exit(False)
        if not end:
            print('Could not find end of dbentry for fragment '\
                '"{}" in line {} of "{}.txt".  Exiting...'.format(fragment, line, db_name))
            sys.exit(False)
        if self.invert:
            atoms = invert_dbatoms_coordinates(atoms)
        return atoms


    def check_consistency(self, fragment_dict, fragment):
        '''
        check if the fragline makes sense and if the fragment_dict is complete
        '''
        for i in fragment_dict:
            if i == 'comment':
                continue
            if not fragment_dict[i]:
                if i == 'head':
                    print('- Restraints in the header of database entry "{}" missing!\n'\
                        ''.format(fragment))
                    return False
                else:
                    print('- Values for "{}" in database entry "{}" missing!'\
                        ''.format(str.upper(i), fragment))
                    return False
        if len(fragment_dict['fragline']) != 8:
            print('- The line starting with "FRAG" in the database entry of {} is '\
                'not correct.\n  Are the cell parameters really correct? '\
                '"FRAG 17 a b c alpha beta gamma"\n'.format(fragment))
            sys.exit(False)
        return True


    def check_db_atom_consistency(self, db, fragment):
        '''This method is for atoms only (without db header)!
          check the db for different types of errors:
        '''
        # check for duplicates:
        atoms = []
        for i in db:
            line = ' '.join(i)
            at = line.split()[0]
            if at in atoms:
                print('duplicate atom {0} in database entry "{1}" '\
                    'found!'.format(at, fragment))
                sys.exit(-1)
            else:
                atoms.append(at)
                continue


    def check_db_header_consistency(self, head, fragment):
        '''
        This method check the db header for consistency
         - before frag only comments or cards are alllowed. in case off error raise warning
           with fragment name and line itself
        :param head: list of strings
        :type head: list
        :param fragment: fragment name
        :type fragment: string
        '''
        for n, line in enumerate(head):
            line = line.upper().split()
            if not line:
                continue
            if line[0] not in SHX_CARDS:  # only the first 4 characters, because SADI_TOL would be bad
                print('Bad line in header of database entry "{}" found!'.format(n, fragment))
                print(' '.join(line))
                sys.exit(False)
        return True



    def get_head_lines(self, fragment, db, line):
        '''
        return the head of the dbentry of the fragment as
        list of strings
        ['RESI CBZ', 'SADI C1 C2 C2 C3 C3 C4 C4 C5 C5 C6 C6 C1',
        'SADI Cl1 C2 Cl1 C6',
        'FLAT Cl1 > C6', 'SIMU Cl1 > C6', 'RIGU Cl1 > C6']
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
            nhead.append(line)
        # nhead is list of strings
        nhead = misc.unwrap_head_lines(nhead)
        self.check_db_header_consistency(nhead, fragment)
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


    def get_atoms_from_fragment(self, fragment):
        '''
        returns the atoms from the dbentry:
        [['O1', '1', '0.01453', '-1.6659', '-0.10966'],
        ['C1', '1', '0.00146', '-0.26814', '-0.06351'], ... ]
        :param fragment: fragment name
        :type fragment: string
        '''
        return self._dbentry_dict[fragment.lower()]['atoms']


    def get_fragline_from_fragment(self, fragment):
        '''
        returns the line with FRAG 17 cell from the dbentry
        '''
        try:
            fragline = self._dbentry_dict[fragment.lower()]['fragline']
        except(KeyError):
            print('Fragment "{}" not found in database!'.format(fragment))
            sys.exit()
        return fragline


    def get_line_number_from_fragment(self, fragment):
        '''
        returns the line lumber from the dbentry
        '''
        return self._dbentry_dict[fragment.lower()]['line']


    def get_head_from_fragment(self, fragment):
        '''
        returns the header of the dbentry of fragment
        '''
        return self._dbentry_dict[fragment.lower()]['head']


    def get_resi_from_fragment(self, fragment):
        '''
        returns the residue name of the dbentry of fragment
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


    def get_db_name_from_fragment(self, fragment):
        '''
        returns the fragment database name of fragment x
        '''
        return self._dbentry_dict[fragment.lower()]['db']


class ImportGRADE():

    def __init__(self, gradefile, invert=False):
        '''
        class to import fragments from GRADE of Global Phasing Ltd.
        :param gradefile:
        :type gradefile:
        :param invert:
        :type invert:
        '''
        self.invert = invert
        self._getdb = ReadDB(dbdir=os.environ["DSR_DB_DIR"],
                 dbnames=["dsr_db.txt", "dsr_user_db.txt"])
        self._db_dir = os.environ["DSR_DB_DIR"]
        self._db_tags = self._getdb.find_db_tags()
        self._gdb = global_DB(invert=False)
        self._db = self._gdb.build_db_dict()
        self._gradefile = gradefile
        gradefiles = self.get_gradefiles()
        self._pdbfile = gradefiles[0]
        self._dfixfile = gradefiles[1]
        self._obpropfile = gradefiles[2]
        self._atoms = self.get_pdbatoms(self._pdbfile)
        self._firstlast = self.get_first_last_atom(self._atoms)
        self._restraints = self.get_restraints()
        self._resi_name = self.get_name_from_obprop(self._obpropfile)
        self._comments = self.get_comments()


    def get_gradefiles(self):
        '''
        returns the .mol2 and .dfx file location
        '''
        grade_base_filename = os.path.splitext(self._gradefile)
        if grade_base_filename[1] == '.tgz':
            gradefile = tarfile.open(self._gradefile)#, encoding="ascii")
        else:
            print('File {} is not a valid file to import from!'.format(grade_base_filename[1]))
            sys.exit(0)
        # names = []
        pdbfile = False
        dfixfile = False
        propfile = False
        for i in gradefile.getnames():
            if i.endswith('.pdb'):
                if re.match(r'.*obabel.*', i):
                    continue
                pdbfile = gradefile.extractfile(i)
            if i.endswith('.dfix'):
                dfixfile = gradefile.extractfile(i)
            if i.endswith('obprop.txt'):
                propfile = gradefile.extractfile(i)
        if not pdbfile:
            print(' .pdb file not found!')
            self.import_error(self._gradefile)
        elif not dfixfile:
            print(' .dfix file not found!')
            self.import_error(self._gradefile)
        elif not propfile:
            print(' obprop.txt file not found!')
            self.import_error(self._gradefile)
        output = []
        output.append(pdbfile.readlines())
        output.append(dfixfile.readlines())
        output.append(propfile.readlines())
        return output


    def get_name_from_obprop(self, obprop):
        '''
        get the fragment name from the obprop.txt file
        :param obprop: file with some information about the molecule
        :type obprop: list of strings
        '''
        regex = re.compile(b'sequence')
        for line in obprop:
            if regex.match(line):
                line = line.split()
                break
            else:
                line = ['found', 'NONE']
        try:
            line[1]
        except(IndexError):
            line = ['found', 'NONE']
        return line[1]


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
        name = 'REM Name: ' + self._resi_name.decode('ascii')
        comments.append(name.split())
        for m in matches:
            for line in self._dfixfile:
                line = line.decode('ascii')
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
        reads the restraints from a dfx-file
        '''
        restraints = []
        for line in self._dfixfile:
            #line = line.decode('ascii')
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
            line = line.decode('ascii')
            line = line.split()
            if line[0] == 'HETATM':
                if line[-1] == 'H':
                    continue
                tmp = []
                tmp.append(line[2])
                tmp.extend(line[6:9])
                atomlines.append(tmp)
        if not atomlines:
            self.import_error()
        if self.invert:
            atomlines = invert_dbatoms_coordinates(atomlines)
        return atomlines


    def bild_grade_db_entry(self):
        '''
        builds a dbentry from the information supplied by GRADEs
        .mol2 and .dfx file
        '''
        db_import_dict = {}
        num = 1
        resi_name = self._resi_name[:3].upper() + str(num)
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
                    print('Importing {} to user database...'.format(name))
                    atomlist = imported_entry[name]['atoms']
                    comment = imported_entry[name]['comment']
                    comment = '\n'.join([' '.join(i) for i in comment if i])
                    head = '\n'.join([' '.join(x) for x in imported_entry[name]['head']])
                    atoms = '\n'.join(['{:<6}  1  {:>8.3f}{:>8.3f}{:>8.3f}'\
                            .format(y[0], float(y[1]), float(y[2]), float(y[3])) for y in atomlist])
                    resi_name = str(name)
                    cell = '  '.join(imported_entry[name]['fragline'])
                    dbentry = '<{}> \n{} \nRESI {} \n{} \n{} \n{} \n</{}>\n''\
                        '.format(resi_name, comment, resi_name, head, cell, atoms, resi_name)
                    f.write(dbentry)
        except(IOError) as e:
            print(e)
            sys.exit(-1)
        # try to write existing dbentrys:
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


    # def read_mol2_file(self, filename):
    #    '''
    #    This methos is deprecated, because we now collect the atoms from
    #    the pdb file!
    #
    #    reads the atom coordiantes from a mol2-file
    #    '''
    #    inputfile = []
    #    print(filename)
    #    try:
    #        gfile = tarfile.open(self._gradefile)
    #        gfile = self._gradefile.extractfile(filename)
    #        inputfile = gfile.readlines()
    #
    #        #with open(filename, 'r') as f:
    #        #    for line in f:
    #        #        inputfile.append(line)
    #    except(IOError) as e:
    #        print(e)
    #        sys.exit(-1)
    #    return inputfile


    def get_name_from_mol2(self):
        '''
        This method is deprecated, because we now collect the name from
        the obprop file!

        get the residue name
        '''
        mol_list = self.read_file(self._molfile)
        regex = r'@<TRIPOS>MOLECULE'
        found = False
        for num, line in enumerate(mol_list):  # @UnusedVariable
            line = line.strip('\n\r')
            if found:
                raw_name = line
                break
            if re.match(regex, line):
                found = True
                continue
        raw_name = raw_name.split('.')[0]
        return raw_name





if __name__ == '__main__':

    gdb = ReadDB()
    dbnames = gdb.find_db_tags()
    # print 'dbcontent:'
    # for i in dbnames:
    #    print ' {:<18}| {:<6}| {:<15}'.format(i[0], i[1], i[2])
    # no valid
    invert = True
    gl = global_DB(invert)
    db = gl.build_db_dict()
    # print db.values()[3]

    # fragment = 'pfanion'
    fragment = 'toluene'
    # fragline = gl.get_fragline_from_fragment(fragment)  # full string of FRAG line
    # dbatoms = gl.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
    dbhead = gl.get_head_from_fragment(fragment)  # this is only executed once
    dbhead = misc.unwrap_head_lines(dbhead)

    # print dbatoms
    # print('residue:', db['toluene']['resi'])
    # print('line of db:', db['toluene']['line'])
    # print('database:', db['toluene']['db'])
    # print(fragline)

    # for i in dbatoms:
    #    print i
    # head = db['toluene']['head']
    print('### DB head:\n')
    for i in dbhead:
        print(i)
    print()
    dbhead = misc.wrap_headlines(dbhead)
    print('### DB head:\n')
    for i in dbhead:
        print(i.strip('\n'))
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
        if i.endswith(('.pdb', '.dfx', '.mol2')):
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

    # for i in mog.get_restraints('./test-data/TOL.dfx'):
    #    print i

    # mog.bild_grade_db_entry()
    # print mog.bild_grade_db_entry('{}.mol2', '{}.dfx').format(dsr_dict[import_grade], dsr_dict[import_grade])
    # mog.write_user_database()
