#/usr/bin/env python
#-*- encoding: utf-8 -*-
#möp

import re
import sys
import fileinput
import itertools
from atoms import Element



SHX_CARDS = ('TITL', 'CELL', 'ZERR', 'LATT', 'SYMM', 'SFAC', 'UNIT', 'LIST', 'L.S.', 'CGLS', 
             'BOND', 'FMAP', 'PLAN', 'TEMP', 'ACTA', 'CONF', 'SIMU', 'RIGU', 'WGHT', 'FVAR', 
             'DELU', 'SAME', 'DISP', 'LAUE', 'REM',  'MORE', 'TIME', 'END',  'HKLF', 'OMIT', 
             'SHEL', 'BASF', 'TWIN', 'EXTI', 'SWAT', 'HOPE', 'MERG', 'SPEC', 'RESI', 'MOVE', 
             'ANIS', 'AFIX', 'HFIX', 'FRAG', 'FEND', 'EXYZ', 'EADP', 'EQIV', 'CONN', 'BIND', 
             'FREE', 'DFIX', 'BUMP', 'SADI', 'CHIV', 'FLAT', 'DEFS', 'ISOR', 'NCSY', 'SUMP',
             'BLOC', 'DAMP', 'STIR', 'MPLA', 'RTAB', 'HTAB', 'SIZE', 'WPDB', 'GRID', 'MOLE', 
             'XNPD', 'REST', 'CHAN', 'FLAP', 'RNUM', 'SOCC', 'PRIG', 'WIGL', 'RANG', 'TANG', 
             'ADDA', 'STAG', 'NEUT', 'ABIN', 'ANSC', 'ANSR', 'NOTR', 'TWST', 'PART') 


def read_resfile(resfile):
    '''read the resfile and return it as list'''
    f = open(resfile)
    reslist = []
    for line in f.readlines():
        reslist.append(line)
    
    f.close()
    return reslist    
    
        

def find_dsr_command(reslist):
    '''find the lines with a DSR command entry and return its line number'''
    regex = 'rem.DSR.put.*'.lower()
    i = []
    print 'searching for following entrys in the resfile:', regex, '\n'
    i = [i for i, l in enumerate(reslist) for m in [re.search(regex, l.lower())] if m]
    try: 
        i[0]
    except(IndexError):
        print 'no DSR command found! Exiting...'
        sys.exit()
    
    if len(i) > 1:
        print 'only one dsr command at once is allowed! Exitig...'
        sys.exit(-1)
    return int(i[0])




def findlinesdb(dbfile):
    '''find the start and endlines of a db entry and return them as list'''
    startline = ''
    endline = ''
    dblines = []
    for dbline in fileinput.input(dbfile):
        #print dbline
        if re.match('<toluene>', dbline):
            startline = fileinput.lineno()
            #print startline, 'start'
        if re.match('</toluene>', dbline):
            endline = fileinput.lineno()
            #print endline, 'end'
        
    dblines.append(startline)
    dblines.append(endline)
    print '\nAn folgenden Positionen wurde was gefunden:', dblines, '\n'''
    return dblines





def extract_db_entry_list(dbfile, dbline):
    '''return a database entry for every dbline given.
    I have to start this funtion for every dsr command once.    
    '''
    entry = []
    startline = dbline[0]
    endline = dbline[1]
    #print startline
    #print endline
    entry.append('\n\n\n rem the following was inserted by DSR: \n')
    db = open(dbfile)
    #put the text between startline and endline into a list:
    for i, line in enumerate(db):
        if i > startline-1:
            entry.append(line)
        if i > endline-3:
            entry.append('  rem end of entry inserted by DSR. \n\n\n')
            break
        
    db.close()
    entry = ''.join(entry)
    return entry




        
def insertentry(dsrline, dbentry, reslist):
    # open the new .res-file and insert the database entry:
    
    newfile = []
    res = []
    nfile = open('p-1-new.res', 'w')

    reslist.insert(dsrline+1,dbentry)  #insert the db entry
    
    ''.join(str(i) for i in reslist)  #convert reslist to string
    
    for item in reslist:
        nfile.write("%s" % item)    #write the new file
    print '\nnew file written!'
    nfile.close()   



    
def find_stringline(reslist, regex):  
    # finds only one occurence:
    '''                  
    for i, s in enumerate(reslist):
        if re.search(regex, s):
            print i
            return i            
    return -1
    '''  
    found = []
    #finds all occourences in the list:  
    for l in (i for i,x in enumerate(reslist) if re.search(regex, x)):
        found.append(l)
    #print found
    return int(found[0])

        


def make_dbentry_with_new_atomtype(dbentry, newlines):
    '''returns the correct dbentry in terms of sfac numbers'''
    nlines = []
    nlines.extend(newlines)
        
    regex = '^([A-Za-z]{1,4}[0-9]{0,3}[A-Za-z]{0,2}\s+[0-9]{1})'
    atomlns = []
    db = dbentry.splitlines()
    
    for l in (i for i,x in enumerate(db) if re.search(regex, x) if not re.search('FRAG', x)):
        atomlns.append(l)
    
    #print atomlns      #lines where atoms are located in the dbentry
    #print db            #splitted dbentry
    #print newlines  #correct atomline in terms of sfac number
    
    for i, x in enumerate(db):
        for y in atomlns:
            if i == y:
                db[i] = nlines[0]
                db[i] = '  '.join(str(i) for i in db[i])
                del nlines[0]
   
    db = '\n'.join(db)
    #print db
    return db 
    
    
    

def get_atomtypes(dbentry, atoms, regex):
    '''find all atoms in a db entry. returns a list of the atomtypes'''
    found = []
    #print atoms
    

    #find lines with atoms and see if they are in the atomlist
    for i, x in enumerate(dbentry.splitlines()):
        z = re.search(regex, x, re.MULTILINE)
        if z:
            l = z.group(1)  # group 1 is the whole atomstring
            y = z.group(2)  # group 2 is only the element
            if y in atoms and len(l) <= 4:
                #print 'found', y, len(l), '\n'
                found.append(y)
    return found
        



def get_sfac_table(reslist):
    '''returns the SFAC table from the reslist'''
    reg = '(SFAC\s+)+([a-zA-Z]{1,2}.*)'
    sfac = ''
    
    for i in reslist:
        m = re.match(reg, i)
        if m:
            sfac = m.group(2)
            
    return sfac.split()
    




def get_db_entry_linenumbers():
    pass
    
    '''das ganze muss so umgeschrieben werden, dass statt des datenbankeintrags 
    ein afix 17 mit den nachfolgenden atomen geschrieben wird.
    Der datenbankeintrag kommt dann oben hinter FVAR'''
    # -create functions for the following:
    # -enumerate atomlines of the dbentry                       ok
    # -enumerate atomnumbers of the DSR-entry                   ok
    # -collect the coordinates of the atoms in the DSR-entry    
    # -generate AFIX with respective atomnumbers
    # -sort the atoms
    # -insert the db entry                                         ok
    # -write a function which reads the scattering factor numbers   ok
    # -find out which scattering factor numbers belong to the       ok
    #  atoms in the database entry an change them accordingly       ok
    #
    # -write a function which is able to generate an AFIX 17 and following atoms.
    #  it should take the number of atoms and their name in the given sort order.
    #
    # -auch ' und " im atomnamen erlauben!
    



def get_atomcoordinates(resfile, in_atoms):
    '''returns the coordiantes of the given atom names as a dictionary'''
    regex = '^([A-Za-z]{1,4}[0-9]{0,3}[A-Za-z]{0,2}\s+[0-9]{1})'
    found_atoms = []
    atom = []
    
    for i in resfile: 
        if re.search(regex, i):            #search atoms
            l = i.strip().split()          #convert to list
            if l[0] not in SHX_CARDS:      #exclude all non-atom cards
                atom.append(l)
    d = {}
    for i in atom:
        x = {i[0]: ' '.join(i[2:5])}
        d.update(x)      #make a dictionary with {'atom name': 'koordinates'}
    
    at = {}
    for x in in_atoms:
        at[x] = d.get(x)  # get each atom coordinate from the d dictionary
    
    return at        #return the in_atoms coordinates as dictionary {atom: 'coordinate'}





def get_db_source_atoms():
    '''returns the source atoms of the db from dsr command'''
    pass



def get_db_target_atoms():
    '''returns the target atoms where the fragment 
    has to be placed from the dsr command'''
    pass






    
def build_afix_entry(resfile,          # results-file
                    db_source_atoms,   # atoms from the db matching the dsr_atoms 
                    res_target_atoms,  # atoms from the dsr command where the fragment should be placed
                    newlines,          # new atoms from the dbentry
                    sfnumbers):        # list of sfac numbers
    '''build an afix entry for the frag fend cards'''
    
    at = []          # list with the dbatom names
    afix_list = []   # final list with atoms, sfac and coordinates
    sfnum = sfnumbers[:]
    
    coord = get_atomcoordinates(resfile, res_target_atoms)
    
    for x in newlines:
        at.append(''.join(str(i) for i in x[0])) #build list with the dbatom names

    for y, i in enumerate(db_source_atoms): 
        if i in at:
            pos = at.index(i) #when atom i of the db_source_atoms is in at, then make pos the index of it
            at[pos] = res_target_atoms[y] #exchange the respective atom in at with the target atom
    
    for n, i in enumerate(at):   # remember: at includes now the target atoms and the dbatoms
        if i in res_target_atoms:
            l = i +'  ' +str(sfnumbers[n]) +'  ' +coord[i]
            afix_list.append(l)               #insert coordinates for the target atoms
        else:
            l =  i +'  ' +str(sfnumbers[n]) +'  0    0    0'  #insert 0 0 0 for all other atoms
            afix_list.append(l)
    
    for i in afix_list:
        print i
    return afix_list
    




def get_sfac_number_of_dbatoms(sfac, types, dbentry):
    '''returns the sfac number of each atom in the db entry as list'''
    #print sfac
    #print types
    sf_numbers = []
    for i in types:   #atomtypes in the db_entry
        for y, x in enumerate(sfac):
            if x == i:
                sf_numbers.append(y+1)
    return sf_numbers
    


def assign_atomtype(db, sfnumbers, types):
    '''returns dbatoms as a list with corrected atomtype'''
    sfnum = sfnumbers[:]
    regex1 = '^([A-Za-z]{1,4}[0-9]{0,3}[A-Za-z]{0,2}\s+[0-9]{1})'
    atomlines = []
    
    for line, i in enumerate(db.splitlines()):  #go through db entry
        try:
            if re.search(regex1, i):               #search atoms
                if not re.search('FRAG', i):       #exclude FRAG line
                    l = i.strip().split()          #convert to list
                    l[1] = sfnum.pop()             #replace scattering factor
                    atomlines.append(l)
        except(IndexError):
            print 'Bad database entry. Exiting...'  #would be nice to have the line number!
            sys.exit(-1)
    return atomlines



if __name__ == '__main__':
    '''main function'''
    
    el = Element()
    
    #start at line beginning, one to 4 characters, zero to 3 digits,
    #one or more whitespaces, one digit, one or more whitespaces, one or more digits,
    #a dot, one or more digits, ... :
    atomregex = '^(([A-Za-z]{1,4})[0-9]{0,3}[A-Za-z]{0,2})\s+[0-9]{1}\s+[0-9]+\.[0-9]+\s+[0-9]+\.[0-9]+\s+[0-9]+\.[0-9]+'
    
    dblines = findlinesdb('database.txt')
    reslist = read_resfile('p-1.res')
    dsrline = find_dsr_command(reslist)
    dbentry = extract_db_entry_list('database.txt', dblines)
    #print dbentry
    
    string = 'FVAR.+[0-9]+'
    fvarline = find_stringline(reslist, string) # line of the FVAR card
    
    
    db_entry_atomtypes = get_atomtypes(dbentry, el.atoms, atomregex)  #atomtypes in the db entry
    sfac_table = get_sfac_table(reslist)   #get the sfac list
    
    sfnumbers = get_sfac_number_of_dbatoms(sfac_table, db_entry_atomtypes, dbentry)
    #print sfnumbers
    newatomlines = assign_atomtype(dbentry, sfnumbers, db_entry_atomtypes)
    #print newatomlines
    
    c_dbentry = make_dbentry_with_new_atomtype(dbentry, newatomlines)
    #print c_dbentry, 'test'
    #insert a frag fend entry in the resfile:
    
    insertentry(fvarline,  #line where the dsr command is found in the resfile
                c_dbentry, #sfac corrected dbentry to insert into resfile
                reslist)   #resfile as list object
    
    #write a function that extracts following lists from the DSR command:
    res_target_atoms = ['C41', 'Q1', 'Q2']
    db_source_atoms = ['C82', 'C83' ,'C84']
    
    #print get_atomcoordinates(reslist, dsr_target_atoms)  #get coordinates of target atoms
    
    
    #build an afix entry which can be inserted into the new res-file:
    build_afix_entry(reslist,           # results-file
                     db_source_atoms,   # atoms from the db matching the dsr_atoms 
                     res_target_atoms,  # atoms from the dsr command where the fragment should be placed in the resfile
                     newatomlines,      # new atoms from the dbentry
                     sfnumbers)
    
    
     
   
