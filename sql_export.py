import sqlite3 as lite

import os

import sys

from atomhandling import get_atomtypes
from constants import SHX_CARDS
from dbfile import ReadDB, global_DB
from atoms import Element

debug = False
#debug = True
if debug:
    import cProfile
    import pstats
    cp = cProfile.Profile()
    cp.enable()

#cp.dump_stats('foo.profile')

gdb = ReadDB()
dbnames = gdb.find_db_tags()
invert = False
gl = global_DB(invert)
db = gl.build_db_dict()
fragment = 'supersilyle'


def get_residue(gl, fragment):
    db_residue_string = gl.get_resi_from_fragment(fragment)
    return db_residue_string


def get_sum_formula(fragment):
    '''
    returns the sum formula of a fragment
    as string like 'SI4 C9 O1'
    :param fragment:
    :type fragment:
    '''
    types = get_atomtypes(db[fragment]['atoms'])
    formula = ''
    for el in set(types):
        num = types.count(el)
        formula+='{}{} '.format(el, num)
    return formula


def get_fragment_atoms_cartesian(fragment):
    '''
    returns the coordinates of the fragment as cartesian coords
    as list of lists [['-2.7538', '15.9724', '22.6810'], ['0.7939', '16.3333', '21.3135'], ...
    :param fragment:
    :type fragment:
    '''
    from export import Export
    gdb = global_DB()
    ex = Export(gdb, invert=False)
    atoms = ex.format_atoms_for_export(fragment)
    coords = []
    for i in atoms:
        co = i.split()[2:5]
        co = ['{:>.4f}'.format(float(x)) for x in co]
        #print(co)
        coords.append(co)
    return coords

cart_coords = get_fragment_atoms_cartesian(fragment)
####################################################
# print(cart_coords)
####################################################

def get_atomic_number(fragment):
    '''
    returns the atomic numbers of the atoms in a fragment in
    same order as the dsr db as list
    :param fragment:
    :type fragment:
    '''
    atnumbers = []
    types = get_atomtypes(db[fragment]['atoms'])
    for i in types:
        el = Element()
        atnumbers.append(el.get_atomic_number(i))
    #print(el.get_element(26))
    return atnumbers

#  numbers = get_atomic_number(fragment)
#############################################
#  print(numbers)
###########################################

def get_atomnames(fragment):
    atoms = db[fragment]['atoms']
    names = []
    for i in atoms:
        names.append(i[0])
    return names

#  atomnames = get_atomnames(fragment)
########################################
#   print(atomnames)
########################################

dbfilename = 'fragment-database.sqlite'
user_dbname = 'user-fragment-database.sqlite'
con = lite.connect(dbfilename)
#conu = lite.connect(user_dbname)
con.text_factory = str
#conu.text_factory = str
cur = con.cursor()
#curu = conu.cursor()


def initialize_db(con, cur):
    print('beginning export')
    con.execute("PRAGMA foreign_keys = ON")
    cur.execute("DROP TABLE IF EXISTS fragment")
    cur.execute("DROP TABLE IF EXISTS atoms")
    cur.execute("DROP TABLE IF EXISTS atom")
    cur.execute("DROP TABLE IF EXISTS FragmentRestraints")
    cur.execute("DROP TABLE IF EXISTS Restraints")
    try:
        cur.execute("DROP INDEX Atoms_FK")
    except:
        pass
    try:
        cur.execute("DROP INDEX Restraint_FK")
    except:
        pass
    try:
        cur.execute("DROP INDEX Fragment_Name")
    except:
        pass
    try:
        cur.execute("DROP INDEX AtomId")
    except:
        pass

    cur.execute('''
                CREATE TABLE Fragment (
                    Id    INTEGER NOT NULL,
                    class  VARCHAR(4),
                    version TEXT,
                    Name    TEXT,
                    Reference    TEXT,
                    comment    TEXT,
                    picture    BLOB,
                    PRIMARY KEY(Id));
                ''')
    cur.execute('''
                CREATE TABLE Atoms (
                    Id    INTEGER NOT NULL,
                    FragmentId    INTEGER NOT NULL,
                    version TEXT,
                    Name    VARCHAR(255),
                    element    VARCHAR(2),
                    x    FLOAT,
                    y    FLOAT,
                    z    FLOAT,
                PRIMARY KEY(Id),
                  FOREIGN KEY(FragmentId)
                    REFERENCES Fragment(Id)
                      ON DELETE CASCADE
                      ON UPDATE NO ACTION);
                ''')
    cur.execute('''
                   CREATE TABLE Restraints (
                      Id INTEGER  NOT NULL,
                      FragmentId  INTEGER NOT NULL,
                      version TEXT,
                      ShelxName CHAR(4),
                      Atoms TEXT,
                    PRIMARY KEY(Id),
                      FOREIGN KEY(FragmentId)
                        REFERENCES Fragment(Id)
                        ON DELETE CASCADE
                        ON UPDATE NO ACTION);
                    ''')

    cur.execute('''
                CREATE INDEX Atoms_FK ON Atoms(FragmentId);
                ''')
    cur.execute('''
                CREATE INDEX Restraint_FK ON Restraints(FragmentId);
                ''')
    cur.execute('''
                CREATE INDEX Fragment_Name ON Atoms(Name);
                ''')
    cur.execute('''
                CREATE INDEX AtomId ON Atoms(Id);
                ''')
    con.execute("PRAGMA foreign_keys = ON")


def fill_fragment_table(table):
    cur.execute('''INSERT INTO fragment (class, version, name, reference, comment, picture)
                                VALUES(?, ?, ?, ?, ?, ?)''', table)

def fill_restraints_table(FragmentId, head):
    for line in head:
        if not line:
            continue
        restr_table = []
        if line[:4] in SHX_CARDS:
            restr_table.append(FragmentId)
            restr_table.append('1')
            restr_table.append(line[:4])
            restr_table.append(line[5:])
        cur.execute('''INSERT INTO Restraints (FragmentId, version, ShelxName, atoms)
                                VALUES(?, ?, ?, ?)''', restr_table)

def fill_atoms(FragmentId, fragment):
    atom_names = get_atomnames(fragment)
    coords = get_fragment_atoms_cartesian(fragment)
    aelements = get_atomic_number(fragment)
    for atname, coord, element in zip(atom_names, coords, aelements):
        x = coord[0]
        y = coord[1]
        z = coord[2]
        table_atoms = (FragmentId, '1', atname, element, x, y, z)
        cur.execute('''INSERT INTO atoms (FragmentId, version, Name, element,
                        x, y, z) VALUES(?, ?, ?, ?, ?, ?, ?)''', table_atoms)

def get_picture(name):
    try:
        fexist = os.stat('./pictures/{}.png'.format(name))
    except OSError:
        fexist = False
    if fexist:
        with open('./pictures/{}.png'.format(name), 'rb') as f:
            return f.read()
    else:
        print('   #####  No picture found!!!')
        sys.exit()

def export_database():
    for fid, fragment in enumerate(db.keys(), 1):
        Name = gl.get_db_name_from_fragment(fragment)
        print(Name)
        #head = '\n'.join(db[fragment]['head'])
        head = db[fragment]['head']
        comment = ' '.join(db[fragment]['comment'])
        #formula = get_sum_formula(fragment)
        reference = gl.get_src_from_fragment(fragment)
        resiclass = get_residue(gl, fragment) 
        picture = get_picture(fragment)
        picture = lite.Binary(picture)
        table_frag = (resiclass, '1', Name, reference, comment, picture)
        #print(table_frag)
        fill_fragment_table(table_frag)
        fill_atoms(fid, fragment)
        fill_restraints_table(fid, head)
    print('\nDatabase exported to "{}"'.format(dbfilename))


    rows = cur.fetchall()
    for row in rows:
        print("{}  {}  {:>8.4f} {:>8.4f} {:>8.4f}".format(*row))
        #print("{}".format(row))


# enable to re-export the db to sql:
initialize_db(con=con, cur=cur)
# Create empty user database:
#initialize_db(con=conu, cur=curu)
export_database()
con.commit()

#t1 = time.clock()
#print_table(con)
#t2 = time.clock()
#print(t2-t1)