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
import sys, re, os
from resfile import ResList
from dsrparse import DSR_Parser
import string
import misc
# all upper case for case insensitivity:
alphabet = [ i for i in string.ascii_uppercase ]
import networkx as nx

# note: parts    1, 2, 3 are _a, _b, _c
# note: residue number 1, 2, 3 are _1, _2, _3
# number first, then part
# part 2 and class 3 = _3b

__metaclass__ = type  # use new-style classes


def remove_duplicate_bonds(bonds):
    '''
    removes duplicates from [(at1, at2, 1.324), (at2, at1, 1.324)]
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


def format_atom_names(atoms, part, resinum):
    '''
    needs a list of atoms ['C1', 'C2', 'O1', ..] witg part number and a residue number
    returns a list with atoms like ['C1_4b', 'C2_4b', 'O1_4b', ..]
    '''
    if resinum:
        pass
    else:
        resinum = ''
    if int(part) > 0:
        partsymbol = alphabet[int(part)-1] # turns part number into a letter
    else:
        partsymbol = ''
    if resinum and partsymbol:
        numpart = '_'+resinum+partsymbol
    if resinum and not partsymbol:
        numpart = '_'+resinum
    else:
        numpart = ''
    # add the _'num''partsymbol' to each atom to be able to find them in the
    # list file:
    atomnames = [i+numpart for i in atoms]
    return atomnames



class ListFile():
    def __init__(self, basefilename):
        self._listfile = basefilename+'.lst'
        self._coord_regex = r'^\s+ATOM\s+x\s+y\s+z'
        self._conntable_regex = r'^.*radii and connectivity'
        self._listfile_list = self.read_lst_file()
        self._single_atom = None
        
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
    
    
    def read_conntable(self):
        '''
        reads the connectivity table from self._listfile_list
        returns a list of all bonded atompairs. Symmetry equivalent atoms 
        are filtered out.
        '''
        symmeq = False
        # find the start of the conntable
        start_line = misc.find_line(self._listfile_list, self._conntable_regex)
        if start_line:
            for num, line in enumerate(self._listfile_list[start_line:]):
                line = line.split()
                # find the end of covalent radii
                if not line:
                    start_line = num+start_line+1
                    break
        connpairs = []
        connlist = []
        for i in self._listfile_list[start_line:]:
            line = i.split()
            if not line:
                break
            if 'found' in line:
                continue
            connlist.append(i.strip(os.linesep).replace('-', '').split())
        for i in connlist:
            atom = i.pop(0)
            for conatom in i:
                if '$' in conatom:
                    symmeq = True
                    continue
                # uppper case for case insensitivity:
                atom = atom.upper()
                conatom = conatom.upper()
                connpairs.append((atom, conatom))
        if symmeq:
            print('\nConnections to symmetry equivalent atoms found.'\
                    ' \nGenerated DFIX restraints might be nonsense!!')
        return connpairs


    def coordinates(self):
        '''
        reads all atom coordinates of the lst-file
        returns a dictionary with {'atom' : ['x', 'y', 'z']}
        '''
        atom_coords = {}
        start_line = int(misc.find_line(self._listfile_list, self._coord_regex))+2
        num = 0
        for line in self._listfile_list[start_line:]:
            line = line.split()
            try:
                line[0]
            except(IndexError):
                break
            xyz = line[1:4]
            xyz = [float(i) for i in xyz]
            atom = {str(line[0]).upper(): xyz}
            atom_coords.update(atom)
        return atom_coords

    
    def get_cell(self):
        '''
        Returns the unit cell parameters from the list file as list:
        ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
        '''
        cell = False
        for line in self._listfile_list:
            if line.startswith(' CELL'):
                cell = line.split()[2:]
                break
        if not cell:
            print('Unable to find unit cell parameters in th res file.')
            sys.exit()
        return cell
    
    @property
    def get_cell_params(self):
        return self.get_cell()
    
    @property
    def get_all_coordinates(self):
        '''
        return all atom coordinates as property 
        {'atom' : ['x', 'y', 'z'], 'atom2' : [...]}
        '''
        return self.coordinates()


    @property
    def get_single_coordinate(self):
        '''
        return the coordinates of a single atom as ['x', 'y', 'z']
        '''
        coord = self.coordinates()
        try:
            return coord[self._single_atom]
        except(KeyError):
            return None

    @get_single_coordinate.setter
    def single_atom(self, atom):
        self._single_atom = atom



class Adjacency_Matrix():
    '''
    returns an adjacence matrix for all atoms in the .lst file.
    edge property is the bond length
    
    needs atoms with numpart
    '''
    
    def __init__(self, atoms, conntable, coords, cell):
        self._atoms = atoms
        self._conntable = conntable
        self._coords = coords
        self._cell = cell
        self.adjmatrix()
        
    def adjmatrix(self):
        '''
        create a distance matrix for the atom coords
        '''
        G=nx.Graph()
        for i in self._conntable:
            atom1 = i[0]
            atom2 = i[1]
            if atom1 in self._atoms:
                coord1 = self._coords[atom1]
                coord2 = self._coords[atom2]
                dist = misc.atomic_distance(coord1, coord2, self._cell)
                G.add_edge(atom1, atom2, weight=dist)
        return G
        
    
    @property
    def get_adjmatrix(self):
        return self.adjmatrix()

    

class Restraints():
    '''
    This class uses the 1,2- and 1,3-bonds from Connections and tries to generate 
    relative distance restraints. 
    Maybe also absolute distance restraints?
    '''
    def __init__(self, coords, G, atoms, cell):
        self.coords = coords
        self.atoms = atoms
        self._cell = cell
        self._G = G


    @property
    def get_12_dfixes(self):
        '''
        returns the requested dfixes als list of strings
        '''
        dfix = []
        for n,i in self._G.adjacency_iter():
            for i, x in list(i.items()):
                dist=x['weight']
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
        '''
        returns the neighbors of the fragment atoms
        '''
        from networkx import exception
        neighbors = []
        nb = None
        for at in atoms:
            try:
                nb = self._G.neighbors(at)
            except(exception.NetworkXError):
                print('Unable to find neighboring atom for "{}"'.format(at))
                #sys.exit()
            neighbors.append([at, nb])
        return(neighbors)
    

    def get_next_neighbors(self):
        '''
        returns the next-neighbors of the fragment atoms
        '''
        nn = []
        nb12 = self.get_neighbors(self.atoms)
        for i in nb12:
            atom1 = i[0]
            bonded = i[1]
            try:
                for n in bonded:
                    pass
            except(KeyError, TypeError):
                print('No atom connections found. DFIX generation not possible!'\
                        ' Check your SHELXL listing file.')
                sys.exit()
            for n in bonded:
                nb = self._G.neighbors(n)
                nb.remove(atom1)
                if not nb:
                    continue
                for at in nb: # nb -> neighbors of n
                    nn.append((atom1, at))
        return(nn)
        
    
    def make_13_dist(self, nn):
        '''
        return 1,3-distance as [('at1', 'at2', 'distance'), ('at1', 'at2', 'distance'), ...]
        needs next neighbors pairs
        '''
        dist_13 = []
        for i in nn:
            atom1 = i[0]
            atom2 = i[1]
            c1 = self.coords[atom1]
            c2 = self.coords[atom2]
            dist_13.append((atom1, atom2, misc.atomic_distance(c1, c2, self._cell)))
        return dist_13


    @property
    def get_formated_13_dfixes(self):
        nextneighbors = remove_duplicate_bonds(self.get_next_neighbors())
        dfixes_13 = self.make_13_dist(nextneighbors)
        dfix_13_format = []
        for i in dfixes_13:
            dfix_13_format.append('DANG {:7}{:7}{:.4f}\n'.format(misc.remove_partsymbol(i[0]), 
                                    misc.remove_partsymbol(i[1]), i[2]))
        return dfix_13_format
        
        



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



if __name__ == '__main__':

    
    from dbfile import global_DB
    from atomhandling import FindAtoms
    from resfile import filename_wo_ending
    import networkx as nx
    from atomhandling import NumberScheme
    res_file = 'p21c.res'
    basefilename = filename_wo_ending(res_file)
    rl = ResList(res_file)
    res_list = rl.get_res_list()
    dsrp = DSR_Parser(res_list, rl)
    dsr_dict = dsrp.parse_dsr_line()
    fragment = 'naphthalene'#dsr_dict['fragment']
    
    gdb = global_DB()
 
    residue = ''
    part = ''
    
    lf = ListFile(basefilename)
    cell = lf.get_cell_params
    lst_file = lf.read_lst_file()
    coords = lf.get_all_coordinates
    conntable = lf.read_conntable()
    dbatoms = gdb.get_atoms_from_fragment(fragment)
    
    num = NumberScheme(res_list, dbatoms, residue)
    fragment_atoms = num.get_fragment_number_scheme()
    #print(numberscheme)
    #fragment_atoms = [i[0] for i in dbatoms]
    fragment_atoms = misc.format_atom_names(fragment_atoms, part, residue)
    print(fragment_atoms, cell)
    am = Adjacency_Matrix(fragment_atoms, conntable, coords, cell)
    G = am.get_adjmatrix
    print(G.nodes())
    print('dihkstra:')
    #print(nx.dijkstra_path(G, 'C1A_B', 'C4A_B'))
    print('\ncycle_basis')
    l = nx.cycle_basis(G)
    # liste der cycles im Graph:
    print(sorted(l))
    re = Restraints(coords, am.get_adjmatrix, fragment_atoms, cell)
    dfixes = re.get_formated_12_dfixes
    dfixes_13 = re.get_formated_13_dfixes
    #print(''.join(dfixes))
    #print(''.join(dfixes_13))
