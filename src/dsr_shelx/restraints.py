# /usr/bin/env python
# -*- encoding: utf-8 -*-
# möp
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <dkratzert@gmx.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#

import os
import re
import string
import sys
from collections import OrderedDict

from atomhandling import get_atomtypes
from elements import get_radius_from_element
from misc import distance, find_line, flatten, get_overlapped_chunks, remove_partsymbol, shift, vol_tetrahedron

# all upper case for case insensitivity:
alphabet = [i for i in string.ascii_uppercase]

# note: parts    1, 2, 3 are _a, _b, _c
# note: residue number 1, 2, 3 are _1, _2, _3
# number first, then part
# part 2 and class 3 = _3b




def remove_duplicate_bonds(bonds):
    """
    removes duplicates from [(at1, at2, 1.324), (at2, at1, 1.324)]
    """
    new_bonds = []
    pairs = []
    for k in bonds:
        pairs.append((k[0], k[1]))
        if (k[1], k[0]) in pairs:
            continue
        new_bonds.append(k)
    return new_bonds


def format_atom_names(atoms, part='', resinum=''):
    """
    needs a list of atoms ['C1', 'C2', 'O1', ..] with part number and a residue number.
    returns a list with atoms like ['C1_4b', 'C2_4b', 'O1_4b', ..]
    :param atoms:  list of plain atom names
    :param part:  string, part number
    :param resinum: string, residue number
    """
    numpart = ''
    part = str(part)
    resinum = str(resinum)
    if not resinum:
        resinum = ''
    try:
        int(part)
        if int(part) > 0:
            partsymbol = alphabet[int(part) - 1]  # turns part number into a letter
        else:
            print('Warning! Part symbol with non-numeric character detected.')
            partsymbol = ''
    except ValueError:
        partsymbol = ''
    if resinum and partsymbol:
        numpart = '_' + resinum + partsymbol
    if not resinum and partsymbol:
        numpart = '_' + partsymbol
    if resinum and not partsymbol:
        numpart = '_' + resinum
    if not resinum and not partsymbol:
        numpart = ''
    # add the _'num''partsymbol' to each atom to be able to find them in the
    # list file:
    atomnames = [i + numpart for i in atoms]
    return atomnames


class Restraints():
    """
    returns an adjacence matrix for all atoms in the .lst file.
    edge property is the bond length

    needs atoms with numpart like C1_2b
    """

    def __init__(self, frag, gdb):
        self.fragment = frag.lower()
        self.gdb = gdb
        self._atoms = self.gdb.get_atoms(fragment=self.fragment, cartesian=False)
        self._atom_names = [i[0] for i in self._atoms]
        self.atom_types = get_atomtypes(self._atoms)
        self.cart_coords = self.gdb.get_coordinates(self.fragment, cartesian=True)
        self._connectivity_table = self.get_conntable_from_atoms(self.cart_coords, self.atom_types, self._atom_names)
        self.coords_dict = self.get_coords_dict()
        self._G = self.get_adjmatrix()

    def get_coords_dict(self):
        """
        Returns an ordered dictionary with coordinates of the fragment
        """
        coords = OrderedDict({})
        for name, co in zip(self._atom_names, self.cart_coords):
            coords[name] = co
        return coords

    def adjmatrix(self):
        """
        create a distance matrix for the atom coordinates
        """
        import networkx as nx
        G = nx.Graph()
        for i in self._connectivity_table:
            atom1 = i[0]
            atom2 = i[1]
            if atom1 in self._atom_names:
                coord1 = self.coords_dict[atom1]
                coord2 = self.coords_dict[atom2]
                dist = distance(coord1[0], coord1[1], coord1[2],
                                coord2[0], coord2[1], coord2[2])
                G.add_edge(atom1, atom2, weight=dist)
        return G

    def get_conntable_from_atoms(self, cart_coords, atom_types, atom_names, extra_param=0.16):
        """
        returns a connectivity table from the atomic coordinates and the covalence
        radii of the atoms.
        TODO:
        - read FREE command from db to control binding here.
        :param cart_coords: cartesian coordinates of the atoms
        :type cart_coords: list
        :param atom_types: Atomic elements
        :type atom_types: list of strings
        :param atom_names: atom name in the file like C1
        :type atom_names: list of strings
        """
        names = []
        for n, i in enumerate(atom_names, 1):
            names.append([n, i])
        conlist = []
        for co1, typ, n1 in zip(cart_coords, atom_types, names):
            for co2, typ2, n2 in zip(cart_coords, atom_types, names):
                if n1 == n2:
                    continue
                ele1_covrad = get_radius_from_element(typ.capitalize())
                ele2_covrad = get_radius_from_element(typ2.capitalize())
                d = distance(co1[0], co1[1], co1[2], co2[0], co2[1], co2[2], round_out=5)
                # print(d, n1, n2, (ele1.covrad+ele2.covrad)+extra_param, '#', ele1.covrad, ele2.covrad)
                # a bond is defined with less than the sum of the covalence
                # radii plus the extra_param:
                if d <= (ele1_covrad + ele2_covrad) + extra_param and d > (ele1_covrad or ele2_covrad):
                    conlist.append([n2[1], n1[1]])
                    if [n1[1], n2[1]] in conlist:
                        continue
        return (conlist)

    def get_adjmatrix(self):
        return self.adjmatrix()

    def get_12_dfixes(self):
        """
        returns the requested dfixes als list of strings
        """
        dfix = []
        for n, i in self._G.adjacency():
            # print(n, i)
            for i, x in list(i.items()):
                dist = x['weight']  # weight of edge is the atomic distance
                atom1 = n
                atom2 = i
                dfix.append((atom1, atom2, dist))
        return dfix

    def get_formated_12_dfixes(self):
        dfix_formated = []
        dfix_restraints = self.get_12_dfixes()
        dfix_restraints = remove_duplicate_bonds(dfix_restraints)
        for n, i in enumerate(dfix_restraints, 1):  # @UnusedVariable
            dfix_formated.append('DFIX {:<7.4f}{:<4s} {:<4s}\n'.format(i[2],
                                                                       remove_partsymbol(i[0]),
                                                                       remove_partsymbol(i[1])))
        return dfix_formated

    def get_neighbors(self, atoms):
        """
        returns the neighbors of the fragment atoms
        """
        from networkx import exception
        neighbors = []
        nb = None
        for at in atoms:
            try:
                nb = self._G.neighbors(at)
            except exception.NetworkXError:
                print('Information: Atom "{}" has no neighbours.'.format(at))
                # sys.exit()
            neighbors.append([at, nb])
        return (neighbors)

    def get_next_neighbors(self):
        """
        returns the next-neighbors of the fragments atoms
        """
        nn = []
        nb12 = self.get_neighbors(self._atom_names)
        for i in nb12:
            atom1 = i[0]
            bonded = list(i[1])
            try:
                for n in bonded:
                    pass
            except(KeyError, TypeError):
                print('No atom connections found. DFIX generation not possible!'
                      ' Check your SHELXL listing file.')
                sys.exit()
            for n in bonded:
                nb = [x for x in self._G.neighbors(n)]
                if not nb:
                    continue
                try:
                    nb.remove(atom1)
                except Exception:
                    # print('Atom {0} has no neighbour.'.format(atom1))
                    pass
                for at in nb:  # nb -> neighbors of n
                    nn.append((atom1, at))
        return nn

    def make_13_dist(self, nn):
        """
        return 1,3-distance as [('at1', 'at2', 'distance'), ('at1', 'at2', 'distance'), ...]
        needs next neighbors pairs
        """
        dist_13 = []
        for i in nn:
            atom1 = i[0]
            atom2 = i[1]
            c1 = self.coords_dict[atom1]
            c2 = self.coords_dict[atom2]
            dist_13.append((atom1, atom2, distance(c1[0], c1[1], c1[2],
                                                   c2[0], c2[1], c2[2])))
        return dist_13

    def make_flat_restraints(self):
        """
        searches for rings in the graph G, splits it in 4-member chunks and tests if
        they are flat: volume of tetrahedron of chunk < 0.1 A-3.
        Additionally, the ring adjacent atoms are added and new chunks created.

        returns list of flat chunks.

        >>> from dbfile import ParseDB
        >>> gdb = ParseDB('./src/dsr_shelx/dsr_db.txt')
        >>> res = Restraints('benzene', gdb)
        >>> sorted(res.make_flat_restraints())
        [['C1', 'C2', 'C3', 'C4'], ['C3', 'C4', 'C5', 'C6']]
        """
        import networkx as nx
        list_of_rings = sorted(nx.cycle_basis(self._G))
        if not list_of_rings:
            return False
        flats = []
        neighbors = []
        newflats = []
        for ring in list_of_rings:
            ring = sorted(ring)
            for atom in ring:
                # lets see if there is a neighboring atom:
                nb = self._G.neighbors(atom)  # [1:]
                for i in nb:
                    if i not in sorted(flatten(list_of_rings)):
                        neighbors.append(i)
            if len(ring) < 4:
                continue  # only proceed if ring is bigger than 3 atoms
            chunks = get_overlapped_chunks(ring, 4)
            for chunk in chunks:
                if self.is_flat(chunk) and chunk not in flats:
                    flats.append(sorted(chunk))
            if not flats:
                return False
            newflats = []
            # check for neighbours and add the to the flat list:
            for chunk in flats:
                newflats.append(chunk)
                for atnum, chunkatom in enumerate(chunk[:]):
                    for nbatom in neighbors:
                        if self.binds_to(nbatom, chunkatom):
                            # add bound atoms near their partners:
                            ch = chunk[:]
                            ch.insert(atnum, nbatom)
                            ch = shift(ch, atnum)
                            H = self._G.subgraph(ch).copy()
                            # Try to delete atoms in the subgraph and test if subgraph divides.
                            # If it not devides, remove the atom unless it is the just added neighbour.
                            for num, i in enumerate(reversed(ch), start=1):
                                # print(ch, nbatom, ch[-num], num, '###')
                                H.remove_node(ch[-num])
                                comp = list(nx.connected_components(H))
                                # check if graph is disconnected now:
                                if len(comp) > 1:
                                    continue
                                else:
                                    # do not delete the just added neighbour:
                                    if ch[-num] in neighbors:
                                        continue
                                    del ch[-num]
                                    break  # finished, go to next flat
                            # only add if it really results in a flat composition:
                            ch.sort()
                            if self.is_flat(ch) and ch not in newflats:
                                newflats.append(ch)
        return newflats

    def binds_to(self, a, b):
        """
        returns True if atom a binds to atoms b
        """
        if [a, b] in self._connectivity_table:
            return True
        else:
            return False

    def is_flat(self, chunk, flat_tresh=0.085):
        # type: (list, float) -> bool
        """
        check if four atoms are flat
        chunk: a list of four atoms coordinates triples
        flat_tresh: Threshold up to whic volume a tetrahedron is flat.
        """
        tetrahedron_atoms = []
        for atom in chunk:
            single_atom_coordinate = self.coords_dict[atom]
            tetrahedron_atoms.append(single_atom_coordinate)
        a, b, c, d = tetrahedron_atoms
        volume = (vol_tetrahedron(a, b, c, d))
        if volume < flat_tresh:
            return True
        else:
            # print('volume of', chunk, 'too big:', volume)
            return False

    def get_formated_flats(self):
        """
        formats the FLAT restraints and removes the part symbol
        """
        flats = self.make_flat_restraints()
        if not flats:
            return []
        flat_format = []
        for i in flats:
            i = [remove_partsymbol(x) for x in i]
            flat_format.append('FLAT {}\n'.format(' '.join(i)))
        return flat_format

    def get_formated_13_dfixes(self):
        nextneighbors = remove_duplicate_bonds(self.get_next_neighbors())
        dfixes_13 = self.make_13_dist(nextneighbors)
        dfix_13_format = []
        for i in dfixes_13:
            dfix_13_format.append('DANG {:<7.4f}{:<4s} {:<4s}\n'.format(i[2],
                                                                        remove_partsymbol(i[0]),
                                                                        remove_partsymbol(i[1])))
        return dfix_13_format


class ListFile():
    def __init__(self, basefilename):
        self._listfile = basefilename + '.lst'
        self._coord_regex = r'^\s+ATOM\s+x\s+y\s+z'
        self._conntable_regex = r'^.*radii and connectivity'
        self._listfile_list = self.read_lst_file()
        self._single_atom = None

    @property
    def listfile_list(self):
        '''
        the current list file as list
        '''
        return self._listfile_list

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

    def read_conntable(self):
        """
        reads the connectivity table from self._listfile_list
        returns a list of all bonded atom pairs. Symmetry equivalent atoms
        are filtered out.
        """
        # find the start of the conntable
        start_line = find_line(self._listfile_list, self._conntable_regex)
        if start_line:
            for num, line in enumerate(self._listfile_list[start_line:]):
                line = line.split()
                # find the end of covalent radii
                if not line:
                    start_line = num + start_line + 1
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
                atom = atom.upper()
                connected_atom = connected_atom.upper()
                connpairs.append((atom, connected_atom))
        return connpairs

    def coordinates(self):
        """
        reads all atom coordinates of the lst-file
        returns a dictionary with {'atom' : ['x', 'y', 'z']}
        """
        atom_coordinates = {}
        start_line = int(find_line(self._listfile_list, self._coord_regex)) + 2
        for line in self._listfile_list[start_line:]:
            line = line.split()
            try:
                line[0]
            except IndexError:
                break
            xyz = line[1:4]
            try:
                xyz = [float(i) for i in xyz]
            except ValueError:
                print('No atoms found in .lst file!')
                sys.exit(0)
            atom = {str(line[0]).upper(): xyz}
            atom_coordinates.update(atom)
        return atom_coordinates

    def get_cell(self):
        """
        Returns the unit cell parameters from the list file as list:
        ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
        """
        cell = False
        for num, line in enumerate(self._listfile_list):
            if line.startswith(' CELL'):
                cell = line.split()[2:]
                try:
                    cell = [float(i) for i in cell]
                except ValueError as e:
                    print('{} \nbad cell parameters in line {} in the list file.'.format(e, num + 1))
                    sys.exit()
                break
        if not cell:
            print('Unable to find unit cell parameters in the .lst file.')
            sys.exit()
        return cell

    @property
    def get_lst_cell_parameters(self):
        """
        Returns the unit cell parameters from the list file as list:
        ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
        """
        return self.get_cell()

    @property
    def get_all_coordinates(self):
        """
        return all atom coordinates as property
        {'atom' : ['x', 'y', 'z'], 'atom2' : [...]}
        """
        return self.coordinates()

    def get_single_coordinate(self, atom):
        """
        return the coordinates of a single atom as ['x', 'y', 'z']
        """
        coord = self.coordinates()
        try:
            return coord[atom]
        except(KeyError):
            return None

    def get_afix_number_of_CF3(self):
        """
        returns the afix number of the atom where SHELXL prints the
        difference density at 15 degree intervals
        """
        regex_atom = r"^\sDifference\selectron\sdensity.*at\s15\sdegree"
        line = find_line(self._listfile_list, regex_atom)
        if not line:
            return False
        at1 = self._listfile_list[line].split('.')[0].split()
        return at1[-1]

    def get_bondvector(self, atom=None):
        """
        get the bond vector in terms of atom names around which SHELXL
        calculates the difference density in 15 degree interval
        Y-Z-F1/F2/F3
        at1 = Y
        at2 = Z
        """
        import networkx as nx
        regex = r'^.*is clockwise looking down'
        line = find_line(self.read_lst_file(), regex)
        if not line and atom:
            conn = self.read_conntable()
            G = nx.Graph(conn)
            nb = G[atom]
            for i in list(nb):
                if i[0] == 'C':
                    return (i, atom)
        if not line:
            return '', ''
        at1 = self._listfile_list[line].split('.')[0].split()[-3]
        at2 = self._listfile_list[line].split('.')[0].split()[-1]
        return at1, at2

    def get_difference_density(self, averaged=False):
        """
        returns the difference density values in 15 degree interval
        if averaged is True, returns the averaged values.
        :param averaged: enable averaged values
        :type averaged: boolean
        """
        regex_atom = r'^.*is clockwise looking down'
        line = find_line(self._listfile_list, regex_atom)
        if not line:
            return False
        if not averaged:
            nums = self._listfile_list[line + 1].split()
        else:
            nums = self._listfile_list[line + 3].split()[4:]
        return [int(i) for i in nums]

    def get_degree_of_highest_peak(self):
        """
        returns the position in degree of the highest peak in the
        difference density. This point is assumed as an atom position.
        """
        dens = self.get_difference_density(averaged=True)
        maximum = dens.index(max(dens)) + 1
        return maximum * 15


class Lst_Deviations():
    """
    reads the deviations of the fitted group from the lst-file
    """

    def __init__(self, lst_file):
        self._lst_file = lst_file
        self._dev = self.find_deviations()

    def find_deviations(self):
        """
        parses the deviations of the fitted group and
        returns a dictionary with the results:
        {'C3': '11.773', 'C2': '7.667', 'C1': '8.761', 'C4': '5.700'}
        """
        regex = re.compile(r'^\s+\*\*\s+Atom.*deviates\sby')
        deviations = {}
        for line in self._lst_file:
            if regex.match(line):
                line = line.split()
                deviations[line[2]] = line[5]
        return deviations

    def print_LS_fit_deviations(self):
        """
        pretty output of the LS-fit deviations
        """
        if self._dev:
            print('\n Fragment fit might have failed!')
            print(' Deviations on fitting group:')
            for i in self._dev:
                print(' {:<4}: {:>5} A'.format(i.strip(' \n\r'), self._dev[i][:4]))


if __name__ == '__main__':
    pass
