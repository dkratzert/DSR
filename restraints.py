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
    needs a list of atoms ['C1', 'C2', 'O1', ..] with part number and a residue number
    returns a list with atoms like ['C1_4b', 'C2_4b', 'O1_4b', ..]
    :param atoms:  list of plain atom names
    :param part:  string, part number
    :param resinum: string, residue number
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
    elif resinum and not partsymbol:
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
            print('Unable to read file {}.'.format(self._listfile))
            sys.exit()
        return listfile


    def read_conntable(self):
        '''
        reads the connectivity table from self._listfile_list
        returns a list of all bonded atom pairs. Symmetry equivalent atoms
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
        connections_list = []
        for i in self._listfile_list[start_line:]:
            line = i.split()
            if not line:
                break
            if 'found' in line:
                continue
            connections_list.append(i.strip(os.linesep).replace('-', '').split())
        for i in connections_list:
            atom = i.pop(0)
            for connected_atom in i:
                if '$' in connected_atom:
                    symmeq = True
                    continue
                # uppper case for case insensitivity:
                atom = atom.upper()
                connected_atom = connected_atom.upper()
                connpairs.append((atom, connected_atom))
        if symmeq:
            print('\nConnections to symmetry equivalent atoms found.'\
                    ' \nGenerated DFIX restraints might be nonsense!!')
        return connpairs


    def coordinates(self):
        '''
        reads all atom coordinates of the lst-file
        returns a dictionary with {'atom' : ['x', 'y', 'z']}
        '''
        atom_coordinates = {}
        start_line = int(misc.find_line(self._listfile_list, self._coord_regex))+2
        for line in self._listfile_list[start_line:]:
            line = line.split()
            try:
                line[0]
            except(IndexError):
                break
            xyz = line[1:4]
            try:
                xyz = [float(i) for i in xyz]
            except(ValueError):
                print('No atoms found in .lst file!')
                sys.exit(0)
            atom = {str(line[0]).upper(): xyz}
            atom_coordinates.update(atom)
        return atom_coordinates


    def get_cell(self):
        '''
        Returns the unit cell parameters from the list file as list:
        ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
        '''
        cell = False
        for num, line in enumerate(self._listfile_list):
            if line.startswith(' CELL'):
                cell = line.split()[2:]
                try:
                    cell = [float(i) for i in cell]
                except(ValueError) as e:
                    print(e, '\nbad cell parameters in line {} in the list file.'.format(num+1))
                    sys.exit()
                break
        if not cell:
            print('Unable to find unit cell parameters in the .lst file.')
            sys.exit()
        return cell

    @property
    def get_lst_cell_parameters(self):
        '''
        Returns the unit cell parameters from the list file as list:
        ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
        '''
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

    needs atoms with numpart like C1_2b
    '''

    def __init__(self, atoms, conntable, coords, cell):
        self._atoms = atoms
        self._connectivity_table = conntable
        self._coordinates = coords
        self._cell = cell
        self.adjmatrix()

    def adjmatrix(self):
        '''
        create a distance matrix for the atom coordinates
        '''
        G=nx.Graph()
        for i in self._connectivity_table:
            atom1 = i[0]
            atom2 = i[1]
            if atom1 in self._atoms:
                coord1 = self._coordinates[atom1]
                coord2 = self._coordinates[atom2]
                dist = misc.atomic_distance(coord1, coord2, self._cell)
                G.add_edge(atom1, atom2, weight=dist)
        return G


    @property
    def get_adjmatrix(self):
        return self.adjmatrix()



class Restraints():
    '''
    This class uses connectivity table from the SHELXL lst file and
    generates 1,2- and 1,3-bond distance restraints.
    '''
    def __init__(self, coordinates, G, atoms, cell):
        '''
        :param coordinates: list of coordinates for all atoms from the lst file
        :param G:  adjacency matrix
        :param atoms:  list, fragment atoms
        :param cell: list, cell parameters
        '''
        self.coordinates = coordinates
        self.atoms = atoms
        cell = [float(i) for i in cell]
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
                dist=x['weight'] # weight of edge is the atomic distance
                atom1 = n
                atom2 = i
                dfix.append((atom1, atom2, dist))
        return dfix


    @property
    def get_formated_12_dfixes(self):
        dfix_formated = []
        dfix_restraints = self.get_12_dfixes
        dfix_restraints = remove_duplicate_bonds(dfix_restraints)
        for n, i in enumerate(dfix_restraints, 1):  # @UnusedVariable
            dfix_formated.append('DFIX {:.4f}  {:7}{:7}\n'.format(i[2], \
                misc.remove_partsymbol(i[0]), misc.remove_partsymbol(i[1])))
        return dfix_formated


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
            c1 = self.coordinates[atom1]
            c2 = self.coordinates[atom2]
            dist_13.append((atom1, atom2, misc.atomic_distance(c1, c2, self._cell)))
        return dist_13


    def get_overlapped_chunks(self, ring, size):
        '''
        returns a list of chunks of size 'size' which overlap with one field.
        If the last chunk is smaller than size, the last 'size' chunks are returned as last chunk.
        '''
        chunks = []
        for i in range(0, len(ring)-size+3, 3):
            chunk = ring[i:i+size]
            if len(chunk) < 4:
                chunk = ring[-size:]
            chunks.append(chunk)
        return chunks


    def make_flat_restraints(self):
        '''
        searches for rings in the graph G, splits it in 4-member chunks and tests if
        they are flat: volume of tetrahedron of chunk < 0.1 A-3.
        returns list of flat chunks.
        '''
        from misc import vol_tetrahedron
        list_of_rings = nx.cycle_basis(self._G)
        #print('The list of rings:', list_of_rings)
        if not list_of_rings:
            return False
        # This creates a list of attached tetrahedron_atoms but is unused atm:
        #attached_atoms = []
        #for ring in list_of_rings:
        #    for atom in ring:
        #        nb = G.neighbors(atom)[1:]
        #        for i in nb:
        #            attached_atoms.append(i)
        #attached_atoms = tuple(set(attached_atoms))
        flats = []
        for ring in list_of_rings:
            if len(ring) < 4:
                continue #wenn ring zu wenig atome hat dann nächsten
            chunks = self.get_overlapped_chunks(ring, 4)
            for p in chunks:
                tetrahedron_atoms = []
                for atom in p:
                    single_atom_coordinate = self.coordinates[atom]
                    tetrahedron_atoms.append(single_atom_coordinate)
                a, b, c, d = tetrahedron_atoms
                volume = (vol_tetrahedron(a, b, c, d, self._cell))
                if volume < 0.1:
                    flats.append(p)
                else:
                    pass
                    #print('volume of', p, 'too big:', volume)
        return flats


    @property
    def get_formated_flats(self):
        '''
        formats the FLAT restraints and removes the part symbol
        '''
        flats = self.make_flat_restraints()
        if not flats:
            return ['']
        flat_format = []
        for i in flats:
            i = [misc.remove_partsymbol(x) for x in i]
            flat_format.append('FLAT {}\n'.format(' '.join(i)))
        return flat_format


    @property
    def get_formated_13_dfixes(self):
        nextneighbors = remove_duplicate_bonds(self.get_next_neighbors())
        dfixes_13 = self.make_13_dist(nextneighbors)
        dfix_13_format = []
        for i in dfixes_13:
            dfix_13_format.append('DANG {:.4f}  {:7}{:7}\n'.format(i[2],
                misc.remove_partsymbol(i[0]), misc.remove_partsymbol(i[1])))
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


    def print_LS_fit_deviations(self):
        '''
        pretty output of the LS-fit deviations
        '''
        if self._dev:
            print('\n Fragment fit might have failed!')
            print(' Deviations on fitting group:')
            for i in self._dev:
                print(' {:<4}: {:>5} A'.format(i.strip(' \n\r'), self._dev[i][:4]))



if __name__ == '__main__':


    from dbfile import global_DB
    #from atomhandling import FindAtoms
    from resfile import filename_wo_ending
    #import networkx as nx
    from atomhandling import NumberScheme
    res_file = 'p21c.res'
    basefilename = filename_wo_ending(res_file)
    rl = ResList(res_file)
    res_list = rl.get_res_list()
    dsrp = DSR_Parser(res_list, rl)
    dsr_dict = dsrp.parse_dsr_line()
    fragment = 'benzene'#dsr_dict['fragment']
    invert = True
    gdb = global_DB(invert)

    residue = '4'
    part = '2'

    lf = ListFile(basefilename)
    cell = lf.get_lst_cell_parameters
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
    print('nodes:', G.nodes())
    # print('dihkstra (kürzester pfad):')
    # print(nx.dijkstra_path(G, 'C1_4C', 'C4_4C'))
    # print('\ncycle_basis')
    l = nx.cycle_basis(G)
    # print('liste der cycles im Graph:')
    print(sorted(l))
    print('end\n')
    restr = Restraints(coords, am.get_adjmatrix, fragment_atoms, cell)
    #dfixes = re.get_formated_12_dfixes
    #dfixes_13 = re.get_formated_13_dfixes
    flats = restr.get_formated_flats
    print(''.join(flats), 'flats')
    #print(''.join(dfixes_13))
