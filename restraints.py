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
import networkx as nx

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
    


class Restraints():
    '''
    This class uses the 1,2- and 1,3-bonds from Connections and tries to generate 
    relative distance restraints. 
    Maybe also absolute distance restraints?
    '''
    def __init__(self):
        pass
        

class Connections():
    '''
    This class parses the shelx .lst file after L.S. 0 refinement and
    puts all 1,2- and 1,3-bonds in a data structure.
    ''' 
    def __init__(self, 
                reslist,   # resfile
                listfile,  # listfile
                dbhead,    # header of dbentry
                atoms,     # atom names for which connections in the list file should be found
                residue):  # residue as number
        self._reslist = reslist
        self._listfile = listfile
        self._dbhead = dbhead
        self._atomnames = atoms[:]
        self._pivot_regex = r'^.*Distance\s+Angles'
    
    
    def get_bond_dist(self):
        '''
        returns a dictionary of connections for each atom in the fragment:
        {'C1_2b': ['C2_2b', '1.3843', 'C6_2b', '1.3849', 'C7_2b', '1.5033', 'C16', '1.9039'], ...
        each value is a set of pairs wit name and distance.
        '''
        # lines is a list where self._pivot_regex is found
        lines = misc.find_multi_lines(self._listfile, self._pivot_regex)
        atomconnections = []
        for num in lines:
            atname = self._listfile[num].split()[0] #atom1
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
                            #angle1 = row[2]
                            bond = ([atname, atname_connect, float(distance)])
                            #bond = ([atname, atname_connect, [float(distance), float(angle1)]])
                            atomconnections.append(tuple(bond))
                        except(IndexError):
                            continue
                    break
        return atomconnections


class Adjacency_Matrix():
    '''
    returns an adjacence matrix for all atoms in the .lst file.
    edge property is the bond length
    
    '''
    def __init__(self, conntable):
        self._conntable = conntable
        
    def build_adjacency_matrix(self):
        '''
        needs pairs of atoms with their distance and residue
        '''
        MG = nx.Graph()
        MG.add_weighted_edges_from(self._conntable)
        return MG

    @property
    def get_adjmatrix(self):
        return self.build_adjacency_matrix()



if __name__ == '__main__':
    from dbfile import global_DB
    from atomhandling import FindAtoms
    options = OptionsParser()
    rl = ResList(options.res_file)
    res_list = rl.get_res_list()
    dsrp = DSR_Parser(res_list, rl)
    dsr_dict = dsrp.parse_dsr_line()
    fragment = dsr_dict['fragment']
    
    gdb = global_DB()
    #dbatoms = gdb.get_atoms_from_fragment(fragment)
    dbhead = gdb.get_head_from_fragment(fragment) 

    lf = ListFile()
    lst_file = lf.read_lst_file()
    #lfd = Lst_Deviations(lst_file)
    #lfd.print_deviations()
        
    fa = FindAtoms(res_list)
    atoms_dict = fa.collect_residues()

# all atom names in the res file:
    res_atoms = misc.get_atoms(res_list)
    res_atoms = [i[0] for i in res_atoms]
   
    con = Connections(res_list, lst_file, dbhead, res_atoms, '2')
    conntable = con.get_bond_dist()
     
    am = Adjacency_Matrix(conntable)
    G = am.get_adjmatrix
    #print(G.nodes())
    print(G.edges(data=True))
    #print(AM.get_adjmatrix())
  #  for n,i in G.adjacency_iter():
  #      for i, x in i.items():
  #          dist=x['weight']
  #          print(dist)
#    for n,nbrs in AM.adjacency_iter():
#        for nbr,eattr in nbrs.items():
#            dist=eattr['weight']
#            #if dist<1.1: print('{}, {}, {:.3f})'.format(n,nbr,dist))
    #print(G.neighbors('Al1'))
    dist='weight'
    #print(G.edge['C1']['C2'][dist])
