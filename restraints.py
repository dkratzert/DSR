
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <daniel.kratzert@ac.uni-freiburg.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
from __future__ import print_function
import sys, re
from resfile import ResList
from options import OptionsParser
from dsrparse import DSR_Parser
from resi import Resi
import string
import misc
from atomhandling import NumberScheme
alphabet = [ i for i in string.ascii_lowercase ]
from dbfile import global_DB


__metaclass__ = type  # use new-style classes

class ListFile():
    def __init__(self):
        options = OptionsParser()
        rl = ResList(options.res_file)
        self._listfile = str(rl.filename_wo_ending(options.res_file))+'.lst'
        
        
    def read_lst_file(self):
        '''
        reads the .lst file and returns it as list. 
        ['Sn1 -       Distance       Angles\n', 'O1_a      2.0997 (0.0088)\n''...']
        
        Connectivity list:
        Sn1 -       Distance       Angles
        O1_a      2.0997 (0.0088) 
        O2        2.1054 (0.0083)   85.87 (0.35)
        O3        2.1176 (0.0082)   83.95 (0.34)  85.10 (0.35)
                    Sn1 -         O1_a          O2 
        
        '''
        listfile = []
        try:
            with open(self._listfile, 'r') as l:
                for line in l:
                    listfile.append(line)
        except(IOError):
            print('Unable to read {}'.format(self._listfile))
            sys.exit()
        return listfile
        
        

class Lst_Deviations():
    '''
    reads the deviations of the fitted group from the lst-file
    '''
    def __init__(self, lst_file):
        self._lst_file = lst_file
        self._dev = self.find_deviations()

    def find_deviations(self):
        '''
        parses the deviations of the fitted group and
        returns a dictionary with the results:
        {'C3': '11.773', 'C2': '7.667', 'C1': '8.761', 'C4': '5.700'}
        '''
        regex = re.compile(r'^\s+\*\*\s+Atom.*deviates\sby')
        deviations = {}
        for line in self._lst_file:
            if regex.match(line):
                line = line.split()
                deviations[line[2]] = line[5]
        return deviations

    
    def print_deviations(self):
        '''
        pretty output of the deviations
        '''
        if self._dev:
            print('\n Fragment fit might have failed.')
            print(' Deviations on fitting group:')
            for i in self._dev:
                print(' {:<4}: {:>5} A'.format(i.strip(' \n\r'), self._dev[i][:4]))
        print('\n')
    

    






# for every atom in fragment:
#   get_neighbours()
#
# note: parts    1, 2, 3 are _a, _b, _c
# note: residue number 1, 2, 3 are _1, _2, _3
# number first, then part
# part 2 and class 3 = _3b
class Connections():
    '''
    This class parses the shelx .lst file after L.S. 0 refinement and
    puts all 1,2- and 1,3-bonds in a data structure.
    
    Simplify __init_!!!!!! This is a mess!!!
    ''' 
    def __init__(self, 
                reslist,   # resfile
                dsr_dict,  # dict of dsr command
                dbhead,    # header of dbentry
                dbatoms,   # final atom names of the dbentry (with new naming sheme)
                fragment,
                residue):  
        self._reslist = reslist
        lf = ListFile()
        self._dbhead = dbhead
        if residue:
            self._resinum = residue
        else:
            self._resinum = ''
        if dsr_dict['part']:
            self._part = dsr_dict['part']
            self._partsymbol = alphabet[int(self._part)-1] # turns part number into a letter
        else:
            self._partsymbol = ''
        if self._resinum and self._partsymbol:
            self._numpart = '_'+self._resinum+self._partsymbol
        else:
            self._numpart = ''
        print('_numpart =',  self._numpart)
##############################################################        
        self._atomnames = dbatoms[:] # instead of dbatoms, use all atoms here to get all connections!
######################################################        
        self._atomnames = [ i+self._numpart for i in self._atomnames] #new names with _numpart notation
        #self._atomnames = ['C1A', 'C2A', 'C3A', 'C4A']
#        self._pivot_regex = r'^\D+.*Distance\s+Angles'
        self._pivot_regex = r'.*Angles.*'
        self._list = lf.read_lst_file()
    
    
    def find_pivotatoms(self):
        '''
        returns a dictionary of connections for each atom in the fragment:
        {'C1_2b': ['C2_2b', '1.3843', 'C6_2b', '1.3849', 'C7_2b', '1.5033', 'C16', '1.9039'], ...
        each value is a set of pairs wit name and distance.
        '''
        #print self._list
        # lines is a list where self._pivot_regex is found
        lines = misc.find_multi_lines(self._list, self._pivot_regex)
        atomconnections = {}
       # print self._atomnames
        #print lines
        for num in lines:
            for atom in self._atomnames:
                #print atom
                if atom in self._list[num]:
                    #print atom
                    atoms = []
                    for i in range(1, 10): # maximum number of connections for each atom
                        try:
                            #print self._list[num+i].split()[1] # 1,2 distance
                            atoms.append(self._list[num+i].split()[0])
                            atoms.append(self._list[num+i].split()[1])
                        except(IndexError):
                            atomconnections[self._list[num].split()[0]] = atoms
                            atoms.pop()
                            atoms.pop()
                            break
        print('test:', atomconnections)
        return atomconnections
    
    
    
    def get_connections(self):
        '''
        get the 1,2 atomic distances of all atoms to all their connections
        '''
        #atoms = []
        atomconnections = {}
        lines = misc.find_multi_lines(self._list, self._pivot_regex)
        for num in lines:
            atoms = []
            atom = self._list[num].split()[0]
            for i in range(1, 10): # maximum number of connections for each atom
                try:
                    #print self._list[num+i].split()[1] # 1,2 distance
                    atoms.append(self._list[num+i].split()[0])
                    atoms.append(self._list[num+i].split()[1])
                except(IndexError):
                    atomconnections[self._list[num].split()[0]] = atoms
                    atoms.pop()
                    atoms.pop()
                    break
        
        print('atomconnections "O4":', atomconnections['O4'])
        
        ##con = Connections(self._reslist, self._dsr_dict, self._dbhead, dbatoms)
        #conntable = self.find_pivotatoms()
        ##print conntable, 'test'
        #moietyconnect = []
        #for i in conntable:
        #    print i,': ' , conntable[i]
        #    for x in conntable[i]:
        #        if x not in conntable:
        #            moietyconnect.append(x)
        #moietyconnect = list(set(moietyconnect)) # these atom(s) connect to the fragment
        ##print moietyconnect, '\nsafsagfesagfr#####################\n'
        ## lets see what else connects to this atom
        #for i in moietyconnect:
        #    if i in conntable:
        #        print i
        #        print conntable
        #        pass

    
    
class Restraints():
    '''
    This class uses the 1,2- and 1,3-bonds from Connections and tries to generate 
    relative distance restraints. 
    Maybe also absolute distance restraints?
    '''
    def __init__(self):
        pass
        



if __name__ == '__main__':
    from dbfile import global_DB
    options = OptionsParser()
    rl = ResList(options.res_file)
    res_list = rl.get_res_list()
    lf = ListFile()
    lst_file = lf.read_lst_file()
    lfd = Lst_Deviations(lst_file)
    lfd.print_deviations()
        
        
        
        
    dsrp = DSR_Parser(res_list)
    dsr_dict = dsrp.parse_dsr_line()
    fragment = dsr_dict['fragment']
    gdb = global_DB()
    dbatoms = gdb.get_atoms_from_fragment(fragment)
    dbhead = gdb.get_head_from_fragment(fragment) 
    num = NumberScheme(res_list, dbatoms, dsr_dict['resi'])
    new_atomnames = num.get_fragment_number_scheme()
    con = Connections(res_list, dsr_dict, dbhead, new_atomnames, 'toluene', '2')
   #    # print 'connections to other atoms:\n'
    pivots = con.find_pivotatoms()
    con.get_connections()
   
    #print 'pivots', pivots
    #here I need the distance of the former atoms list and make them in an adjacence matrix
    
    
    #for i in pivots:
    #    print i,':' , ' '.join(con.find_pivotatoms()[i])
    #   
    #for values in pivots.keys():
    #    for i in pivots.get(values):
    #        if i in pivots.keys() and i not in pivots.keys():
    #            print i, 'bindet an', values
