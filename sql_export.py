
"""
This file is used to export the DSR database to sqlite for Olex2.
Every fragemnt needs a picture with the same name as the fragment tag in ./pictures
"""
import sqlite3 as lite
import os
import sys

import misc
from constants import SHX_CARDS
from dbfile import global_DB

misc.remove_file('./fragment-database.sqlite')
gl = global_DB(invert=False, maindb="./olex_dsr_db.txt", userdb='no userdatabase!')
db = gl.build_db_dict()


dbfilename = 'fragment-database.sqlite'
con = lite.connect(dbfilename)
con.text_factory = str
cur = con.cursor()


def get_fragment_atoms_cartesian(fragment):
    """
    returns the coordinates of the fragment as cartesian coords
    as list of lists [['-2.7538', '15.9724', '22.6810'], ['0.7939', '16.3333', '21.3135'], ...
    :param fragment:
    :type fragment:
    """
    from export import Export
    ex = Export(gl, invert=False)
    atoms = ex.format_atoms_for_export(fragment)
    coords = []
    for i in atoms:
        co = i.split()[2:5]
        co = ['{:>.4f}'.format(float(x)) for x in co]
        # print(co)
        coords.append(co)
    return coords


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


def fill_restraints_table(fragmentid, head):
    for line in head:
        if not line:
            continue
        restr_table = []
        if line[:4] in SHX_CARDS:
            restr_table.append(fragmentid)
            restr_table.append('1')
            restr_table.append(line[:4])
            restr_table.append(line[5:])
        cur.execute('''INSERT INTO Restraints (fragmentid, version, ShelxName, atoms)
                                VALUES(?, ?, ?, ?)''', restr_table)


def fill_atoms(FragmentId, fragment):
    atom_names = gl.get_atomnames(fragment)
    coords = get_fragment_atoms_cartesian(fragment)
    aelements = gl.get_atomic_numbers(fragment)
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
        print('   #####  No picture  for {} found!!!'.format(name))
        sys.exit()


def export_database():
    for fid, fragment in enumerate(db.keys(), 1):
        Name = gl.get_name_from_fragment(fragment)
        print("Exporting {}: {}: {}".format(fid, fragment, Name))
        # head = '\n'.join(db[fragment]['head'])
        head = db[fragment]['head']
        comment = ' '.join(db[fragment]['comment'])
        # formula = db.get_sum_formula(fragment)
        reference = gl.get_src_from_fragment(fragment)
        resiclass = gl.get_resi_from_fragment(fragment)
        picture = get_picture(fragment)
        picture = lite.Binary(picture)
        table_frag = (resiclass, '1', Name, reference, comment, picture)
        fill_fragment_table(table_frag)
        fill_atoms(fid, fragment)
        fill_restraints_table(fid, head)
    print('\nDatabase exported to "{}"'.format(dbfilename))

    rows = cur.fetchall()
    for row in rows:
        print("{}  {}  {:>8.4f} {:>8.4f} {:>8.4f}".format(*row))
        # print("{}".format(row))



# enable to re-export the db to sql:
initialize_db(con=con, cur=cur)
# Create empty user database:
# initialize_db(con=conu, cur=curu)
export_database()
con.commit()

# t1 = time.clock()
# print_table(con)
# t2 = time.clock()
# print(t2-t1)
