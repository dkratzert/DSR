#-*- encoding: utf-8 -*-
#m�p
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <daniel.kratzert@ac.uni-freiburg.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
from __future__ import print_function
import string
from constants import *
import re
import os

alphabet = string.ascii_uppercase

__metaclass__ = type  # use new-style classes

def get_atoms(atlist):
    '''returns all atoms found in the input as list'''
    atoms = []
    for i in atlist:
        if re.search(atomregex, str(i)):        # search atoms
            l = i.split()[:5]              # convert to list and use only first 5 columns
            if l[0].upper() not in SHX_CARDS:      # exclude all non-atom cards
                atoms.append(l)
    return atoms


def ll_to_string(inputlist):
    '''converts list of list to string with four whitespaces between each list element'''
    inputlist = [list(map(str, i)) for i in inputlist]
    newlist = []
    for i in inputlist:
        newlist.append('   '.join(i).rstrip())
    string = '\n'.join(newlist)
    #string = '\n'.join('\t'.join(map(str, l)).rstrip() for l in inputlist).expandtabs(6).rstrip()
    return string


def multiline_test(line):
    '''
    test if the current line is a multiline with "=" at the end
    '''
    line = line.rpartition('=')  # partition the line in before, sep, after
    line = ''.join(line[0:2])    # use all including separator
    line = line.rstrip()         # strip spaces
    if line.endswith('='):
        return True
    else:
        return False


def find_line(inputlist, regex):
    '''
    returns the index number of the line where regex is found in the __inputlist
    '''
    for i, string in enumerate(inputlist):
        if re.match(regex, string, re.IGNORECASE):
            return i      # returns the index number if regex found
    return False          # returns False if no regex found (xt solution has no fvar)


def find_multi_lines(inputlist, regex):
    '''
    returns the index number of all lines where regex is found in the inputlist
    '''
    reg = re.compile(regex)
    foundlist = []
    for i, string in enumerate(inputlist):
        if reg.match(string):
            foundlist.append(i)      # returns the index number if regex found
        else:
            continue
    return foundlist


def remove_file(filename, exit_dsr=False, terminate=False):
    '''
    removes the file from disk
    program exits when exit is true
    platon gets terminated if terminate is true
    '''
    if os.path.isfile(filename):
        try:
            os.remove(filename)
        except(WindowsError, OSError):
            #print 'unable to cleanup ins {} files!'.format(file)
            if terminate:
                pgrogname=terminate
                pgrogname.terminate()
            if exit_dsr:
                sys.exit(0)
    

def wrap_headlines(dbhead):
    import textwrap
    for num, line in enumerate(dbhead):
        line = textwrap.wrap(line, 77, subsequent_indent = '  ')
        if len(line) > 1:
            newline = []
            for n, l in enumerate(line):
                if n < len(line)-1:
                    l = l+' =\n'
                newline.append(l)
            dbhead[num] = ' '.join(newline)
    for num, line in enumerate(dbhead):
        line = ' '.join(line.strip().split(' '))
        dbhead[num] = line+'\n'
    return dbhead



def unwrap_head_lines(headlines):
    '''
    if a line is wrapped like "SADI C1 C2 =\n  C3 C4" or "SADI C1 C2=\n  C3 C4"
    this function returns "SADI C1 C2 C3 C4"
    '''
    import constants
    tmp = ''
    new_head = []
    for line in headlines:
        line = line.strip(' \n=')+' '
        tmp = tmp+line
    line = tmp.split()
    for n, i in enumerate(line):
        if i in constants.SHX_CARDS:
            line[n] = '\n'+line[n]
    new_head = ' '.join(line).strip().split('\n')
    return new_head


# this is deprecated:
def get_replace_mode(dsr_string):
    '''returns the put/replace keyword'''
    dsr = makelist(dsr_string)
    action = dsr[2]
    if action == 'REPLACE':
        return True
    else:
        return False


def makelist(string):
    '''returns an upper-case list from a text string'''
    stringlist = []
    for i in string.split():
        stringlist.append(i.upper())
    return stringlist

    
def which(name, flags=os.X_OK):
    '''Search PATH for executable files with the given name.
    
    On newer versions of MS-Windows, the PATHEXT environment variable will be
    set to the list of file extensions for files considered executable. This
    will normally include things like ".EXE". This fuction will also find files
    with the given name ending with any of these extensions.

    On MS-Windows the only flag that has any meaning is os.F_OK. Any other
    flags will be ignored.
    '''
    result = []
    #exts = filter(None, os.environ.get('PATHEXT', '').split(os.pathsep))
    exts = ['.exe', '.EXE']
    path = os.environ.get('PATH', None)
    if path is None:
        return []
    for p in os.environ.get('PATH', '').split(os.pathsep):
        p = os.path.join(p, name)
        if os.access(p, flags):
            result.append(p)
        for e in exts:
            pext = p + e
            if os.access(pext, flags):
                result.append(pext)
    return result    
    
    
def zero(m,n):
    '''
    Create zero matrix of dimension m,n
    '''
    new_matrix = [[0 for row in range(n)] for col in range(m)]
    return new_matrix

 
def matrix_mult(matrix1,matrix2):
    '''
    Multiplies matrix1 with matrix2.
    Independent from numpy, but slow.
    '''
    if len(matrix1[0]) != len(matrix2):
        # Check matrix dimensions
        print('Matrices must be m*n and n*p to multiply!')
    else:
        # Multiply if correct dimensions
        new_matrix = zero(len(matrix1),len(matrix2[0]))
        for i in range(len(matrix1)):
            for j in range(len(matrix2[0])):
                for k in range(len(matrix2)):
                    new_matrix[i][j] += matrix1[i][k]*matrix2[k][j]
    return new_matrix





def format_atom_names(atoms, part, resinum):
    '''
    needs a list of atoms ['C1', 'C2', 'O1', ..] witg part number and a residue number
    returns a list with atoms like ['C1_4b', 'C2_4b', 'O1_4b', ..]
    '''
    if not resinum:
        resinum = ''
    try:
        int(part)
        partsymbol = alphabet[int(part)-1] # turns part number into a letter
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


def remove_partsymbol(atom):
    '''
    strips the part symbol like C1_4b from an atom name
    '''
    if '_' in atom:
        prefix = atom.split('_')[0]
        suffix = atom.split('_')[-1].strip(string.ascii_letters)
        if not suffix:
            atom = prefix
        else:
            atom = prefix+'_'+suffix
    else:
        atom = atom+'_0'
    return atom


def atomic_distance(p1, p2, cell):
    '''
    p1 and p2 are x, y , z coordinates as list ['x', 'y', 'z']
    cell are the cell parameters as list: ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
    returns the distance between the two points.
    ''' 
    from math import cos, sin, sqrt, radians
    cell = [float(y) for y in cell]
    a , b, c =  cell[:3]
    al = radians(cell[3])
    be = radians(cell[4])
    ga = radians(cell[5])
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    dx = (x1-x2)
    dy = (y1-y2)
    dz = (z1-z2)
    dsq = (a*dx)**2+\
          (b*dy)**2+\
          (c*dz)**2+\
          2*b*c*cos(al)*dy*dz+\
          2*dx*dz*a*c*cos(be)+\
          2*dx*dy*a*b*cos(ga)
    return(sqrt(dsq))


def frac_to_cart(frac_coord, cell):
    '''
    Converts fractional coordinates to cartesian coodinates
    '''
    from math import cos, sin, sqrt, radians
    a, b, c, alpha, beta, gamma = cell
    x, y, z = frac_coord
    alpha = radians(alpha)
    beta  = radians(beta)
    gamma = radians(gamma)
    cosastar = (cos(beta)*cos(gamma)-cos(alpha))/(sin(beta)*sin(gamma))
    sinastar = sqrt(1-cosastar**2)
    Xc = a*x + (b*cos(gamma))*y + (c*cos(beta))*z
    Yc = 0   + (b*sin(gamma))*y + (-c*sin(beta)*cosastar)*z
    Zc = 0   +  0               + (c*sin(beta)*sinastar)*z
    return (Xc, Yc, Zc)


def determinante(a):
    '''
    return determinant of 3x3 matrix
    '''
    return (a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2])
           -a[1][0] * (a[0][1] * a[2][2] - a[2][1] * a[0][2])
           +a[2][0] * (a[0][1] * a[1][2] - a[1][1] * a[0][2]))


def subtract(a, b):
    '''
    subtract vector b from vector a
    '''
    return (a[0] - b[0],
            a[1] - b[1],
            a[2] - b[2])


def vol_tetrahedron(a, b, c, d, cell):
    '''
    returns the volume of a terahedron spanned by four points:
    e.g. A = (3, 2, 1), B = (1, 2, 4), C = (4, 0, 3), D = (1, 1, 7)

            |u1 u2 u3|
    v = 1/6*|v1 v2 v3|
            |w1 w2 w3|
    
    AB = (1-3, 2-2, 4-1) = (-2, 0, 3)
    AC = ...
    AD = ...
    
    V = 1/6[u,v,w]
    
              |-2,  0, 3|
    [u,v,w] = | 1, -2, 2| = 24-3-12 = 5
              |-2, -1, 6|
        
    V = 1/6*5
    '''
    a = [float(i) for i in a]
    b = [float(i) for i in b]
    c = [float(i) for i in c]
    d = [float(i) for i in d]
    A = frac_to_cart(a, cell)
    B = frac_to_cart(b, cell)
    C = frac_to_cart(c, cell)
    D = frac_to_cart(d, cell)
    
    AB = subtract(A, B)
    AC = subtract(A, C)
    AD = subtract(A, D)
    D = determinante([AB, AC, AD])
    volume = abs((D/6))

    return volume
    


#def bond_angle(dx1, dx2, cell):
#    pass
#    # cos(phi) = {a**2*dxr*dxs + b**2*dyr*dys + c**2 * dzr*dzs + b*c*cos(alpha)(dyr*dzs + dys*dzr) +  
#   #      c*a*cos(beta)(dzr*dxs + dzs*xr) + a*b*cos(gamma)*(dxr*dys + dxs*dyr) } / r s
#    vector1 = (2,3,5)
#    vector2 = (3,4,6)
#    def dot():
#        sum(p*q for p,q in zip(vector1, vector2))
#    def normalize(v):
#        vmag = magnitude(v)
#        return [ v[i]/vmag  for i in range(len(v)) ]
#    def magnitude(v):
#        return math.sqrt(sum(v[i]*v[i] for i in range(len(v))))
#    def get_angle(self, list):
#        """Get angle formed by three atoms.
#
#        calculate angle between the vectors list[1]->list[0] and
#        list[1]->list[2], where list contains the atomic indexes in
#        question."""
#        # normalized vector 1->0, 1->2:
#        v10 = self.positions[list[0]] - self.positions[list[1]]
#        v12 = self.positions[list[2]] - self.positions[list[1]]
#        v10 /= np.linalg.norm(v10)
#        v12 /= np.linalg.norm(v12)
#        angle = np.vdot(v10, v12)
#        angle = np.arccos(angle)
#        return angle

if __name__ == '__main__':
    from resfile import ResList, ResListEdit
    from atomhandling import FindAtoms
    from dsrparse import DSR_Parser
    import math as m
    import sys
    from dbfile import global_DB
    res_file = 'p21c.res'
    res_list = ResList(res_file)
    reslist = res_list.get_res_list()
    find_atoms = FindAtoms(reslist)
    rle = ResListEdit(reslist, find_atoms)
    dsrp = DSR_Parser(reslist, rle)
    
    gdb = global_DB(self.invert)
    db = gdb.build_db_dict()
    fragment = 'PFAnion'
    fragline = gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
    #dbatoms = gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
    dbhead = gdb.get_head_from_fragment(fragment)        # this is only executed once
    resi = True #gdb.get_resi_from_fragment(fragment)
    uhead = unwrap_head_lines(dbhead)
    print(uhead)
    sys.exit()
    
    
    
    
    cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
    # CF3:
    a = (0.281319, 0.368769, 0.575106)
    b = (0.352077, 0.314955, 0.582945)
    c = (0.191896, 0.365437, 0.617058)
    d = (0.358359, 0.417667, 0.594043)
    print('volume of CF3-group:')
    print(vol_tetrahedron(a, b, c, d, cell))
    # Benzene:
    # orig: a = (0.838817,   0.474526,   0.190081)
    a = (0.838817,   0.484526,   0.190081) # a ist um 0.01 ausgelenkt
    b = (0.875251,   0.478410,   0.256955)
    c = (0.789290,   0.456520,   0.301616)
    d = (0.674054,   0.430194,   0.280727)
    print('volume of Benzene ring atoms:')
    print(vol_tetrahedron(a, b, c, d, cell))    
    
    head = ['FLAT C C1 C10 C11 C12 C13 C2 C3', 'FLAT C C1 C10 C11 C12 C13 C2 C3 C4 C5 C6 C7 C8 C9 CL C4 C5 C6 C7 C8 C9 CL C4 C5 C6 C7 C8 C9 CL C4 C5 C6 C7 C8 C9 CL C8 C9 CL C4 C5 C6 C7 C8 C9 CL C4 C5 C6 C7 C8 C9 CL  C8 C9 CL C4 C5 C6 C7 C8 C9 CL C4 C5 C6 C7 C8 C9 CL C8 C9 CL C4 C5 C6 C7 C8 C9 CL C4 C5 C6 C7 C8 C9 CLx']
    
    whead = wrap_headlines(head)
    #print(whead)
    uhead = unwrap_head_lines(whead)
    #print(uhead)
    sys.exit()
    dsr_string = dsrp.find_dsr_command(line=True).lower()

    regex = 'Q2.*'
    print('found regex in line', find_line(reslist, regex))
    regex2 = '^H[0-9]+\s+'
    multi = find_multi_lines(reslist, regex2)
    print('found multiline in lines', multi)
    if not multi:
        print('nix gefunbden!!!!!!!!')
    print('replace mode:', get_replace_mode(dsr_string))
    print('#'+reslist[123].strip('\n')+'#')
    print('multiline?', multiline_test(reslist[123]))
    print('#'+reslist[24].strip('\n')+'#')
    print('multiline?', multiline_test(reslist[24]))
    
    
    cell90 = (1, 1, 1, 90, 90, 90)
    cell90 = [ float(i) for i in cell90 ]
    coord1 = (-0.186843,   0.282708,   0.526803) # C5
    #                                              -2.741   5.912  10.774
    coord2 = (-0.155278,   0.264593,   0.600644) # C7
    # 1.573 A                                      -2.520   5.533  12.289
    
    N1 = frac_to_cart(coord1, cell)
    N2 = frac_to_cart(coord2, cell)
    print(N1, '-2.741   5.912  10.774 must be same')
    print(N2, '-2.520   5.533  12.289')
    x1 = float(N1[0])
    y1 = float(N1[1])
    z1 = float(N1[2])
    
    x2 = float(N2[0])
    y2 = float(N2[1])
    z2 = float(N2[2])
    
    d = m.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    
    print('\ndist_frac_to_cart           : {:.3f}'.format(d))
    print('dist_frac_to_cart_atomicdist: {:.3f}'.format(atomic_distance(N1, N2, cell90)))
    print('atomicdist_direct:          : {:.3f}'.format(atomic_distance(coord1, coord2, cell)))
    print('korrekte dist:              : 1.573 A\n')
    
    coo = frac_to_cart((0.312, 0.37, 0.754), (10.5086, 20.9035, 20.5072, 90, 94.13, 90))
    #print coo
    print('{:8.6} {:8.6} {:8.6}'.format(*coo), 'neue koordianten')
