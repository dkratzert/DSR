#/usr/bin/env python
#-*- encoding: utf-8 -*-
#möp
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
#import networkx as nx

# note: parts    1, 2, 3 are _a, _b, _c
# note: residue number 1, 2, 3 are _1, _2, _3
# number first, then part
# part 2 and class 3 = _3b

__metaclass__ = type  # use new-style classes

class ListFile():
    def __init__(self):
        options = OptionsParser()
        rl = ResList(options.res_file)
        self._listfile = str(rl.filename_wo_ending(options.res_file))+'.lst'
        self._coord_regex = r'^\s+ATOM\s+x\s+y\s+z'
        self._listfile_list = self.read_lst_file()
        
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
    
    
    def coordinates(self, listfile):
        '''
        reads the atom coordinates of the lst-file
        '''
        atom_coords = {}
        start_line = int(misc.find_line(listfile, self._coord_regex))+2
        #print(listfile)
        num = 0
        for line in listfile[start_line:]:
            line = line.split()
            try:
                line[0]
            except(IndexError):
                break
            atom = {str(line[0]): line[1:4]}
            atom_coords.update(atom)
        return atom_coords
        
    @property
    def get_coordinates(self):
        '''
        return the coordinates as property
        '''
        return self.coordinates(self._listfile_list )
        
    
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
    


class Restraints():
    '''
    This class uses the 1,2- and 1,3-bonds from Connections and tries to generate 
    relative distance restraints. 
    Maybe also absolute distance restraints?
    '''
    def __init__(self, G):
        self._G = G
    
    @property
    def get_12_dfixes(self):
        '''
        returns the requested dfixes als list of strings
        '''
        single_dfix = []
        all_dfixes = []
        for n,i in self._G.adjacency_iter():
            for i, x in i.items():
                dist=x['1,2-dist']
                single_dfix = ' '.join([n, i, dist])
                all_dfixes.append(single_dfix+'\n')
                print(single_dfix)
            return all_dfixes
    
    

class Connections():
    '''
    This class parses the shelx .lst file after L.S. 0 refinement and
    puts all 1,2-bonds in a data structure.
    ''' 
    def __init__(self, 
                reslist,   # resfile
                listfile,  # listfile
                dbhead,    # header of dbentry
                atoms,     # atom names for which connections in the list file should be found
                part,      # part number
                residue):  # residue as number
        self._reslist = reslist
        self._listfile = listfile
        self._dbhead = dbhead
        self._pivot_regex = r'^.*Distance\s+Angles'
        if residue:
            self._resinum = residue
        else:
            self._resinum = ''
        if part:
            self._part = part
            self._partsymbol = alphabet[int(self._part)-1] # turns part number into a letter
        else:
            self._partsymbol = ''
        if self._resinum and self._partsymbol:
            self._numpart = '_'+self._resinum+self._partsymbol
        else:
            self._numpart = ''
        #print('_numpart =',  self._numpart)
        self._atomnames = [i+self._numpart for i in atoms]
        self._atomnames_resi = [i+self._resinum for i in atoms]
       # print(self._atomnames)
    
    @property
    def get_numpart(self):
        return self._numpart
    
    @property
    def get_atoms_plusresi(self):
        return self._atomnames_resi
    
    def get_bond_dist(self):
        '''
        returns a list of connections for each atom in the fragment:
        ['C1, 'C2' {'1,2-dist': distance}, 'O1, 'Al1' {'1,2-dist': distance}]
        each value is a set of pairs wit name and distance.
        '''
        # lines is a list where self._pivot_regex is found
        lines = misc.find_multi_lines(self._listfile, self._pivot_regex)
        atomconnections = []
        for num in lines:
            atname = self._listfile[num].split()[0] #atom1
            #print(self._listfile[num].split())
            for atom in self._atomnames:
                if atom in self._listfile[num]:
                    bond = []
                    for i in range(1, 10): # maximum number of connections for each atom
                        row = self._listfile[num+i].split()
                        try:
                            atname_connect = row[0] #bonded atom
                            distance = row[1]
                            if distance == '-':
                                break
                            bond = [atname, atname_connect, {'1,2-dist':float(distance)}]
                            atomconnections.append(tuple(bond))
                        except(IndexError):
                            continue
                    break
        return atomconnections



class Nextneighbors():
    '''
    Nextneighbors for 1,3-distances
    '''
    
    
    def get_13_dist():
        '''
        
        '''
        def __init__(self):
            pass
    

class Adjacency_Matrix():
    '''
    returns an adjacence matrix for all atoms in the .lst file.
    edge property is the bond length
    
    [('O1_4b', 'C1_4b', {'1,2-dist': 1.3986}), ... ]
    '''
    
    def __init__(self, conntable, residue):
        self._conntable = conntable
        self._residue = residue
        
    def build_adjacency_matrix(self):
        '''
        needs pairs of atoms with their distance and residue
        '''
        import networkx as nx
        MG = nx.Graph()
        MG.add_edges_from(self._conntable)
        return MG

    @property
    def get_residue(self):
        return self._residue
    
    @property
    def get_adjmatrix(self):
        return self.build_adjacency_matrix()




if __name__ == '__main__':
    # -list of bonds from Covalent radii and connectivity table for fragment
    # -coordinates for each atom (node) in above list
    # -adjacency matrix for each bond
    # -calculate weight (distcance) for each bond (edge)
    # -then 1,3 distance -> next neighbors and weight again
    
    from dbfile import global_DB
    from atomhandling import FindAtoms
    options = OptionsParser()
    rl = ResList(options.res_file)
    res_list = rl.get_res_list()
    dsrp = DSR_Parser(res_list, rl)
    dsr_dict = dsrp.parse_dsr_line()
    fragment = 'oc(cf3)3'#dsr_dict['fragment']
    
    gdb = global_DB()
    #dbatoms = gdb.get_atoms_from_fragment(fragment)
    dbhead = gdb.get_head_from_fragment(fragment) 

    lf = ListFile()
    lst_file = lf.read_lst_file()
    
    #lfd = Lst_Deviations(lst_file)
    #lfd.print_deviations()
        
    fa = FindAtoms(res_list)
    atoms_dict = fa.collect_residues()
    dbatoms = gdb.get_atoms_from_fragment(fragment)
    dbatoms = [i[0] for i in dbatoms]
    
    coords = lf.get_coordinates
    #print(coords['Al1'])
    
# all atom names in the res file:
    #res_atoms = misc.get_atoms(res_list)
    #res_atoms_list = [i[0] for i in res_atoms]
    
    con = Connections(res_list, lst_file, dbhead, dbatoms, dsr_dict['part'], '3')
    conntable = con.get_bond_dist()
    #print(conntable)
    am = Adjacency_Matrix(conntable, '3')
    G = am.get_adjmatrix
    
    #print(G.nodes())
    #print(G.edges(data=True))
    #print(atoms_dict)
 #   for i in (dbatoms):
 #       #print(i[0])
  #      print(G.edges(i+'_4b', data=True))
    #print(AM.get_adjmatrix())
    print()
    dfix = []
    for n,i in G.adjacency_iter():
        #print(i)
        for i, x in i.items():
            dist=x['1,2-dist']
            dfix.append([n, i, dist])
    #print(dfix)
    
    
    def get_neighbors(atoms):
        neighbors = []
        for p in atoms:
            #print(p[0])#, p[1])
            #next = (G.neighbors(p[0]), 'test1', p[0]) #p[0] macht 1,3 zu test2-p[0]
            #print(G.neighbors(p[1]), 'test2\n')
            nb = G.neighbors(p[1])
            try:
                nb.remove(p[0])
            except(ValueError):
                pass
            except(AttributeError):
                pass
    #        print(p[0], nb) #1, 3
            if not nb:
                pass
            else:
                neighbors.append([p[0], nb])
            #print(neighbors)
        return(neighbors)
    
    print('rtest')
    
    nb = get_neighbors(dfix)
    print(nb)
    
    
    length = len(nb)
    for n, item in enumerate(nb):
        piv = item[0]
        if n == length/2:
            break
        for k in item[1]:
            pass
            #print(piv, k)
            #print(piv, k, fa.get_atomcoordinates([piv]))
        

    
    #misc.at_distance(misc.cell(res_list))
    
    #jetzt matrix für 1,3 machen und für jedes atom die koordinaten
    
    
    # nb enthält jetzt jeweils ein atom pro listenlement und seine dazugehörigen 1,3-stehenden
    
    
    #for i in nb:
    #    print(G.neighbors(i))
    
    
#    for n,nbrs in AM.adjacency_iter():
#        for nbr,eattr in nbrs.items():
#            dist=eattr['weight']
#            #if dist<1.1: print('{}, {}, {:.3f})'.format(n,nbr,dist))
    #print(G.neighbors('Al1'))
    #dist='weight'
    #print(G.edge['C1']['C2'][dist])
