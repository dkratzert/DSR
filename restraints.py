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
from dsrparse import DSR_Parser
import string
import misc
from collections import OrderedDict
from atomhandling import get_atomtypes
from misc import distance, vol_tetrahedron
from elements import ELEMENTS
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


def format_atom_names(atoms, part='', resinum=''):
    '''
    needs a list of atoms ['C1', 'C2', 'O1', ..] with part number and a residue number.
    returns a list with atoms like ['C1_4b', 'C2_4b', 'O1_4b', ..]
    :param atoms:  list of plain atom names
    :param part:  string, part number
    :param resinum: string, residue number
    '''
    part = str(part)
    resinum = str(resinum)
    if not resinum:
        resinum = ''
    try:
        int(part)
        if int(part) > 0:
            partsymbol = alphabet[int(part)-1] # turns part number into a letter
        else:
            print('Warning! Part symbol with non-numeric character detected.')
            partsymbol = ''
    except(ValueError):
        partsymbol = ''
    if resinum and partsymbol:
        numpart = '_'+resinum+partsymbol
    if not resinum and partsymbol:
        numpart = '_'+partsymbol
    if resinum and not partsymbol:
        numpart = '_'+resinum
    if not resinum and not partsymbol:
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



class Restraints():
    '''
    returns an adjacence matrix for all atoms in the .lst file.
    edge property is the bond length

    needs atoms with numpart like C1_2b
    '''

    def __init__(self, fragment, gdb):
        fragment = fragment.lower()
        self.gdb = gdb
        self.db = self.gdb.build_db_dict()
        self._atoms = [i[0] for i in self.db[fragment]['atoms']]
        self._cell = self.gdb.get_unit_cell(fragment)
        self.fragment = fragment
        cart_coords = self.get_fragment_atoms_cartesian()
        self.cart_coords = [[float(y) for y in i] for i in cart_coords]
        self.atom_types = get_atomtypes(self.db[fragment]['atoms'])
        self._connectivity_table = self.get_conntable_from_fragment()
        self.coords_dict = self.get_coords_dict()
        self._G = self.get_adjmatrix()
        
    def get_coords_dict(self):
        coords = OrderedDict({})
        for name, co in zip(self._atoms, self.cart_coords):
            coords[name] = co
        return coords

    def get_fragment_atoms_cartesian(self):
        '''
        returns the coordinates of the fragment as cartesian coords
        as list of lists [['-2.7538', '15.9724', '22.6810'], ['0.7939', '16.3333', '21.3135'], ...
        :param fragment:
        :type fragment:
        '''
        from export import Export
        ex = Export(self.fragment, self.gdb, False)
        atoms = ex.format_atoms_for_export()
        coords = []
        for i in atoms:
            coords.append(i.split()[2:5])
        return coords

    def adjmatrix(self):
        '''
        create a distance matrix for the atom coordinates
        '''
        G=nx.Graph()
        #print(self._atoms, self.coords_dict)
        for i in self._connectivity_table:
            atom1 = i[0]
            atom2 = i[1]
            if atom1 in self._atoms:
                coord1 = self.coords_dict[atom1]
                coord2 = self.coords_dict[atom2]
                dist = misc.distance(coord1[0], coord1[1], coord1[2], \
                                     coord2[0], coord2[1], coord2[2])
                G.add_edge(atom1, atom2, weight=dist)
        return G

    def get_conntable_from_fragment(self):
        '''
        returns a connectivity table from the atomic coordinates and the covalence
        radii of the fragment atoms.
        '''
        names = []
        for n, i in enumerate(self._atoms, 1):
            names.append([n, i])
        conlist = []
        for co1, typ, n1 in zip(self.cart_coords, self.atom_types, names):
            for co2, typ2, n2 in zip(self.cart_coords, self.atom_types, names):
                ele1 = ELEMENTS[typ.capitalize()]
                ele2 = ELEMENTS[typ2.capitalize()]
                d = distance(co1[0], co1[1], co1[2], co2[0], co2[1], co2[2], round_out=5)
                if d <= (ele1.covrad+ele2.covrad)+0.05 and d != 0.0:
                    if n1 == n2:
                        continue
                    conlist.append([n2, n1])
                    if [n1, n2] in conlist:
                        continue
                    #print('{}--{}: {}'.format(n1, n2, d))
        conlist = [(i[0][1], i[1][1]) for i in conlist]
        #print('####', conlist)
        return (conlist) 

    def get_adjmatrix(self):
        return self.adjmatrix()


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

    def get_formated_12_dfixes(self):
        dfix_formated = []
        dfix_restraints = self.get_12_dfixes()
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
        nb12 = self.get_neighbors(self._atoms)
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
            c1 = self.coords_dict[atom1]
            c2 = self.coords_dict[atom2]
            dist_13.append((atom1, atom2, distance(c1[0], c1[1], c1[2], 
                                                   c2[0], c2[1], c2[2])))
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
        
        first add neighbor atoms to neighbors
        check if original rings are flat, if flat check if ring with neighbor
        is flat, if yes, add this chunk minus first atom
        '''
        list_of_rings = nx.cycle_basis(self._G)
        #print('The list of rings:', list_of_rings)
        if not list_of_rings:
            return False
        flats = []
        neighbors = []
        for ring in list_of_rings:
            for atom in ring:
                # lets see if there is a neighboring atom:
                nb = self._G.neighbors(atom)[1:]
                for i in nb:
                    if not i in ring:
                        neighbors.append(i)
            if len(ring) < 4:
                continue #wenn ring zu wenig atome hat dann nächsten
            chunks = self.get_overlapped_chunks(ring, 4)
            for chunk in chunks:
                if self.is_flat(chunk):
                    flats.append(chunk[:])
                for i in neighbors:
                    chunk.append(i)
                    del chunk[0]
                    del neighbors[0]
                    if self.is_flat(chunk):
                        flats.append(chunk)
        return flats


    def is_flat(self, chunk):
        '''
        check if four atoms are flat
        '''
        tetrahedron_atoms = []
        for atom in chunk:
            single_atom_coordinate = self.coords_dict[atom]
            tetrahedron_atoms.append(single_atom_coordinate)
        a, b, c, d = tetrahedron_atoms
        volume = (vol_tetrahedron(a, b, c, d))
        if volume < 0.1:
            return True
        else:
            print('volume of', chunk, 'too big:', volume)
            return False
            
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
    fragment = 'OC(CF3)3'
    fragment= fragment.lower()
    invert = True
    gdb = global_DB(invert)

    residue = '4'
    part = '2'

    dbatoms = gdb.get_atoms_from_fragment(fragment)

    num = NumberScheme(res_list, dbatoms, residue)
    fragment_atoms = num.get_fragment_number_scheme()
    #print(numberscheme)
    #fragment_atoms = [i[0] for i in dbatoms]
    fragment_atoms = misc.format_atom_names(fragment_atoms, part, residue)
    #print(fragment_atoms)

    restr = Restraints(fragment, gdb)
    dfix_12 = restr.get_formated_12_dfixes()
    dfix_13 = restr.get_formated_13_dfixes()
    for i in dfix_12:
        print(i.strip('\n'))
    for i in dfix_13:
        print(i.strip('\n'))
    flats = restr.get_formated_flats()
    print('flats:\n'+''.join(flats))
    #print(''.join(dfixes_13))

    def make_eadp(self):
        '''
        find atoms nearby and make them eadp
        '''
        for at in self.all_atoms:
            pass