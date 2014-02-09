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
        #print('\n')
    


class Restraints():
    '''
    This class uses the 1,2- and 1,3-bonds from Connections and tries to generate 
    relative distance restraints. 
    Maybe also absolute distance restraints?
    '''
    def __init__(self, conntable, residue, res_list, fa):
        from resfile import get_cell
        self.fa = fa
        self.cell = get_cell(res_list)
        am = Adjacency_Matrix(conntable, residue)
        self._G = am.get_adjmatrix
        self.nb = self.get_neighbors(self.get_12_dfixes)


    @property
    def get_12_dfixes(self):
        '''
        returns the requested dfixes als list of strings
        '''
        dfix = []
        for n,i in self._G.adjacency_iter():
            for i, x in list(i.items()):
                dist=x['1,2-dist']
                atom1 = n
                atom2 = i
                dfix.append((atom1, atom2, dist))
        return dfix
    
    
    @property
    def get_formated_12_dfixes(self):
        dfix_format = []
        dfix = self.get_12_dfixes
        dfix = remove_duplicate_bonds(dfix)
        for n, i in enumerate(dfix, 1):
            dfix_format.append('DFIX {:7}{:7}{:.4f}\n'.format(misc.remove_partsymbol(i[0]), 
                                misc.remove_partsymbol(i[1]), i[2]))
        return dfix_format


    def get_neighbors(self, atoms):
        neighbors = []
        for p in atoms:
            nb = self._G.neighbors(p[1])
            try:
                nb.remove(p[0])
            except(ValueError):
                pass
            except(AttributeError):
                pass
            if not nb:
                pass
            else:
                atom1 = p[0]
                atom2 = nb
                #print(p)
                neighbors.append([atom1, atom2])
                #for i in nb:
                #    print(nx.shortest_path(G, source=atom1,target=i))
            #print(neighbors)
        return(neighbors)


    def get_coordpairs(self, nb):
        '''
        coordinate pairs for 1,3-distances
        '''
        coordpairs = []
        for n, item in enumerate(nb):
            piv = item[0]
            for k in item[1]:
                coordpairs.append((self.fa.get_atomcoordinates([k]), self.fa.get_atomcoordinates([piv])))
        return coordpairs
    
    
    def make_13_dist(self, coordpairs):
        '''
        return 1,3-distance as [('at1', 'at2', 'distance'), ('at1', 'at2', 'distance'), ...]
        '''
        distpairs_13 = []
        for at1, at2 in coordpairs:
            for coord1, coord2 in zip(list(at1.keys()), list(at2.keys())):
                c1 = [ float(i) for i in at1[coord1] ]
                c2 = [ float(i) for i in at2[coord2] ]
                atom1 = misc.remove_partsymbol(list(at1.keys())[0])
                atom2 = misc.remove_partsymbol(list(at2.keys())[0])
                distpairs_13.append((atom1, atom2, misc.at_distance(c1, c2, self.cell)))
        distpairs_13 = remove_duplicate_bonds(distpairs_13)
        return distpairs_13


    def get_13_dist(self):
        coordpairs = self.get_coordpairs(self.nb)
        try:
            dist_13 = self.make_13_dist(coordpairs)
        except(AttributeError):
            print('No DFIX restraints could be inserted.')
            sys.exit()
        return dist_13

    @property
    def get_formated_13_dfixes(self):
        dfixes_13 = self.get_13_dist()
        dfix_13_format = []
        for n, i in enumerate(dfixes_13, 1):
            dfix_13_format.append('DANG {:7}{:7}{:.4f}\n'.format(misc.remove_partsymbol(i[0]), 
                                    misc.remove_partsymbol(i[1]), i[2]))
        return dfix_13_format
        
        

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
                part,      # part number as string
                residue):  # residue number as string
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
        self.atoms = [i[0] for i in atoms]
        self._atomnames = [i+self._numpart for i in self.atoms]
        self._atomnames_resi = [i+'_'+str(self._resinum) for i in self.atoms]

    
    @property
    def get_numpart(self):
        return self._numpart
    
    @property
    def get_atoms_plusresi(self):
        return self._atomnames_resi
    
    def get_bond_dists(self):
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
    
    




def remove_duplicate_bonds(bonds):
    '''
    removes duplicates from at1 at2 1.324, at2 at1 1.324
    '''
    #keys = set()
    new_bonds = [] #(dict([(pair[0], pair) for pair in reversed(bonds)]).values())
    pairs = []
    for k in bonds:
        pairs.append((k[0], k[1]))
        if (k[1], k[0]) in pairs:
            continue
        new_bonds.append(k)
    return new_bonds




        
    
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
        MG = nx.Graph()
        MG.add_edges_from(self._conntable)
        return MG

    @property
    def get_residue(self):
        return self._residue
    
    @property
    def get_adjmatrix(self):
        return self.build_adjacency_matrix()

    @property
    def get_graph(self):
        return self.build_adjacency_matrix


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
    fa = FindAtoms(res_list)
    atoms_dict = fa.collect_residues()
    dbatoms = gdb.get_atoms_from_fragment(fragment)
    coords = lf.get_coordinates
    residue = '4'
    con = Connections(res_list, lst_file, dbhead, dbatoms, '2', residue)
    conntable = con.get_bond_dist()
    re = Restraints(conntable, residue, res_list)
    dfixes = re.get_formated_12_dfixes
    dfixes_13 = re.get_formated_13_dfixes
    print(''.join(dfixes))
    print(''.join(dfixes_13))
