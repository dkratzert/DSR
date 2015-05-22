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
import re
import os
from constants import atomregex, SHX_CARDS
import math
from math import cos, sqrt, radians, sin
import shutil

alphabet = string.ascii_uppercase

__metaclass__ = type  # use new-style classes

reportlog = 'dsr_bug_report.log'

def checkFileExist(filename):
    '''
    Check if a file exists and has some content.
    A file zize above 0 returns True.
    A non-existing file returns False.
    A file size of 0 retrurns 'zero'.
    '''
    filesize = False
    status = False
    try:
        filesize = int(os.stat(str(filename)).st_size)
    except:
        print('File "{}" not found!'.format(filename))
        status = False
    if isinstance(filesize, int) and filesize > 0:
        status = True
    if isinstance(filesize, int) and filesize == 0:
        status = 'zero'
    return status


def get_atoms(atlist):
    '''
    returns all atoms found in the input as list of lists
    input:
    ['F8    4    0.349210   0.073474   0.519443 -21.00000   0.03106', '...']
    output:
    [['O1', '3', '-0.01453', '1.66590', '0.10966'], ['C1', '1', '-0.00146', '0.26814', '0.06351'], ...
    '''
    atoms = []
    try:
        atlist[0]
    except:
        return []
    if isinstance(atlist[0], list):
        for i in atlist:
            atoms.append(i[:5])
        return atoms
    for i in atlist:
        if re.match(atomregex, str(i)):        # search atoms
            l = i.split()[:5]              # convert to list and use only first 5 columns
            if l[0].upper() not in SHX_CARDS:      # exclude all non-atom cards
                atoms.append(l)
    return atoms


def flatten(nested):
    '''
    flattens a nested list
    '''
    result = []
    try:
        # dont iterate over string-like objects:
        try: nested + ''
        except(TypeError): pass
        else: raise TypeError
        for sublist in nested:
            for element in flatten(sublist):
                result.append(element)
    except(TypeError):
        result.append(nested)
    return result

def sortedlistdir(directory):
    '''
    returns a sorted list of files in directory directory.
    :param directory: directory
    :type directory: string
    :param cmpfunc: compare funtion to sort
    :type cmpfunc: string
    '''
    dirlist = os.listdir(directory)
    dirlist.sort()
    return dirlist


def find_line_of_residue(reslist, resinumber):
    '''
    Returns the line number where residue n appears in the reslist.
    :param reslist: res file as list like:
                    ['C1 1 -0.00146 0.26814 0.06351 11.00 0.05',
                    'RESI 4 BENZ',
                    'C2 1 -1.13341 -0.23247 -0.90730 11.00 0.05',]
    :type reslist: list of lists
    :param resinumber: residue number
    :type resinumber: string
    :return n: integer
    :return line: string
    '''
    for n, line in enumerate(reslist):
        if line.upper().startswith('RESI'):
            if line.split()[1] == str(resinumber):
                return [n, line]


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
    :param line: 'O1 3 -0.01453 1.66590 0.10966 11.00 0.05 ='
    :type line: string
    '''
    line = line.rpartition('=')  # partition the line in before, sep, after
    line = ''.join(line[0:2])    # use all including separator
    line = line.rstrip()         # strip spaces
    if line.endswith('='):
        return True
    else:
        return False


def find_line(inputlist, regex, start=None):
    '''
    returns the index number of the line where regex is found in the inputlist
    if stop is true, stop searching with first line found
    :param inputlist: list of strings
    :type inputlist: list
    :param regex: regular expression to search
    :type regex: string
    :param start: line number where to start the search
    '''
    if start:
        start = int(start)
        inputlist_slice = inputlist[start:]
        for i, string in enumerate(inputlist_slice, start):
            if re.match(regex, string, re.IGNORECASE):
                return i      # returns the index number if regex found
    else:
        for i, string in enumerate(inputlist):
            if re.match(regex, string, re.IGNORECASE):
                return i      # returns the index number if regex found
    return False          # returns False if no regex found (xt solution has no fvar)


def find_multi_lines(inputlist, regex):
    '''
    returns the index number of all lines where regex is found in the inputlist
    ! this method is case insensitive !
    '''
    reg = re.compile(regex, re.IGNORECASE)
    foundlist = []
    for i, string in enumerate(inputlist):
        if reg.match(string):
            foundlist.append(i)      # returns the index number if regex found
        else:
            continue
    return foundlist


def remove_file(filename, exit_dsr=False, terminate=False):
    '''
    removes the file "filename" from disk
    program exits when exit is true
    platon gets terminated if terminate is true
    '''
    if os.path.isfile(filename):
        try:
            os.remove(filename)
        except(WindowsError, OSError):  # @UndefinedVariable
            print('can not delete {}'.format(file))
            #print 'unable to cleanup ins {} files!'.format(file)
            if terminate:
                pgrogname=terminate
                pgrogname.terminate()
            if exit_dsr:
                sys.exit(0)

def copy_file(source, target):
    '''
    Copy a file from source to target. Source can be a single file or 
    a directory. Target can be a single file or a directory. 
    :param source: list or string
    :param target: string
    TODO: implement list as source
    '''
    #target_file = os.path.basename(target)
    target_path = os.path.dirname(target)
    source_file = os.path.basename(source)
    if isinstance(source, (list, tuple)):
        print('can not copy a list.') 
    if not os.path.exists(target_path):
        try:
            os.makedirs(target_path)
        except(IOError, OSError):
            print('Unable to create directory {}.'.format(target_path))
    try:
        shutil.copyfile(source, target)
    except(IOError):
        print('Unable to copy {}.'.format(source_file))


def make_directory(dirpath):
    '''
    create a directory with all subdirs from the last existing path
    :param dirpath: string
    '''
    try:
        os.makedirs(dirpath)
    except(IOError, OSError):
        print('Unable to create directory {}.'.format(dirpath))
        #sys.exit(False)


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
    if a line is wrapped like "SADI C1 C2 =\n", "  C3 C4" or "SADI C1 C2=\n", "  C3 C4"
    this function returns "SADI C1 C2 C3 C4"
    '''
    import constants
    tmp = ''
    # it is faster with this loop:
    eq = False
    for line in headlines:
        if '=' in line:
            eq = True
            break
    if not eq:
        return headlines
    for line in headlines:
        line = line.strip(' \n=')
        line = line.replace('=', ' ')
        tmp = tmp+' '+line
    line = tmp.split()
    for n, i in enumerate(line):
        if i[:4] in constants.SHX_CARDS:
            line[n] = '\n'+line[n]
    new_head = ' '.join(line).strip().split('\n')
    new_head = [i.rstrip(' ') for i in new_head]
    return new_head


def makelist(string):
    '''
    returns an upper-case list from a text string
    '''
    stringlist = []
    for i in string.split():
        stringlist.append(i.upper())
    return stringlist


def which(name, flags=os.X_OK):
    '''
    Search PATH for executable files with the given name.

    On MS-Windows the only flag that has any meaning is os.F_OK. Any other
    flags will be ignored.
    '''
    result = []
    #exts = filter(None, os.environ.get('PATHEXT', '').split(os.pathsep))
    exts = ['.exe', '.EXE']
    path = os.getenv('PATH', None)
    if path is None:
        return []
    for p in os.getenv('PATH', '').split(os.pathsep):
        p = os.path.join(p, name)
        if os.access(p, flags):
            result.append(p)
        for e in exts:
            pext = p + e
            if os.access(pext, flags):
                result.append(pext)
    return result


def remove_partsymbol(atom):
    '''
    strips the part symbol like C1_4b from an atom name
    :param atom: 'C1_4b'
    :type atom: string
    '''
    if '_' in atom:
        prefix = atom.split('_')[0]
        suffix = atom.split('_')[-1].strip(string.ascii_letters)
        if not suffix:
            atom = prefix
        else:
            atom = prefix+'_'+suffix
    #else:
    #    atom = atom+'_0  ! This restraint might be inaccurate.'
    return atom


def atomic_distance(p1, p2, cell):
    '''
    p1 and p2 are x, y , z coordinates as list ['x', 'y', 'z']
    cell are the cell parameters as list: ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
    returns the distance between the two points.
    '''
    from math import cos, sqrt, radians
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
    :param frac_coord: [float, float, float]
    :param cell:       [float, float, float, float, float, float]
    '''
    #from math import cos, sin, sqrt, radians
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
    return (round(Xc, 8), round(Yc, 8), round(Zc, 8))


def cart_to_frac(cart_coord, cell):
    '''
    converts cartesian coordinates to fractional coordinates
    :param frac_coord: [float, float, float]
    :param cell:       [float, float, float, float, float, float]
    '''
    a, b, c, alpha, beta, gamma = cell
    X, Y, Z = cart_coord
    alpha = radians(alpha)
    beta  = radians(beta)
    gamma = radians(gamma)
    cosastar = (cos(beta)*cos(gamma)-cos(alpha))/(sin(beta)*sin(gamma))
    sinastar = sqrt(1-cosastar**2)
    z = Z/(c*sin(beta)*sinastar) 
    y = (Y-(-c*sin(beta)*cosastar)*z)/(b*sin(gamma))
    x = (X-(b*cos(gamma))*y-(c*cos(beta))*z)/a
    return (round(x, 8), round(y, 8), round(z, 8))
    
def zero(m,n):
    '''
    Create zero matrix of dimension m,n
    :param m: integer
    :param n: integer
    '''
    new_matrix = [[0 for row in range(n)] for col in range(m)]  # @UnusedVariable
    return new_matrix


def matrix_mult(matrix1,matrix2):
    '''
    Multiplies matrix1 with matrix2.
    Independent from numpy, but slow.
    [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]
    '''
    if len(matrix1[0]) != len(matrix2):
        # Check matrix dimensions
        print('Matrices must be m*n and n*p to multiply!')
        return False
    else:
        # Multiply if correct dimensions
        new_matrix = zero(len(matrix1),len(matrix2[0]))
        for i in range(len(matrix1)):
            for j in range(len(matrix2[0])):
                for k in range(len(matrix2)):
                    new_matrix[i][j] += matrix1[i][k]*matrix2[k][j]
    return new_matrix

def mm(A, B):
    '''
    matrix multiplication (might be faster than above)
    '''
    return [[sum(x * B[i][col] for i,x in enumerate(row)) for col in range(len(B[0]))] for row in A]


def matrix_mult_vector(A, v):
    '''
    multiplies a matrix with a vector
    | 00 01 02 | |0|
    | 10 11 12 |*|1|
    | 20 21 22 | |2|
    '''
    if len(A[0]) != len(v):
        print('Size of matrix and vector not equal.')
        raise Exception
    #a = A[0][0]*v[0]+A[0][1]*v[1]+A[0][2]*v[2]
    #b = A[1][0]*v[0]+A[1][1]*v[1]+A[1][2]*v[2]
    #c = A[2][0]*v[0]+A[2][1]*v[1]+A[2][2]*v[2]
    new_matrix = zero(len(A[0]),0)
    for i in range(len(A)):
        for k in range(len(v)):
            new_matrix[i] += A[i][k]*v[k]
    return new_matrix

def translate_coords(coords):
    '''
    translate coordinates
    '''
    pass


def determinante(a):
    '''
    return determinant of 3x3 matrix
    '''
    return (a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2])
           -a[1][0] * (a[0][1] * a[2][2] - a[2][1] * a[0][2])
           +a[2][0] * (a[0][1] * a[1][2] - a[1][1] * a[0][2]))


def cross_vec(a, b):
    '''
    Cross product of two 3D vectors
    '''
    assert len(a) == len(b) == 3, 'For 3D vectors only'
    a1, a2, a3 = a
    b1, b2, b3 = b
    return (a2*b3 - a3*b2, a3*b1 - a1*b3, a1*b2 - a2*b1)


def subtract_vect(a, b):
    '''
    subtract vector b from vector a
    :param a: [float, float, float]
    :param b: [float, float, float]
    '''
    return (a[0] - b[0],
            a[1] - b[1],
            a[2] - b[2])


def transpose(a):
    '''
    transposes a matrix
    '''
    return zip(*a)

def norm_vec(a):
    '''
    returns a normalized vector
    '''
    l = sqrt(a[0]**2 + a[1]**2 + a[2]**2)
    return (a[0]/l, a[1]/l, a[2]/l)
    

def vol_tetrahedron(a, b, c, d, cell=None):
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
    A = [float(i) for i in a]
    B = [float(i) for i in b]
    C = [float(i) for i in c]
    D = [float(i) for i in d]
    if cell:
        A = frac_to_cart(a, cell)
        B = frac_to_cart(b, cell)
        C = frac_to_cart(c, cell)
        D = frac_to_cart(d, cell)
    AB = subtract_vect(A, B)
    AC = subtract_vect(A, C)
    AD = subtract_vect(A, D)
    D = determinante([AB, AC, AD])
    print(D, '###')
    volume = abs((D/6))
    return volume

def vol_unitcell(a, b, c, al, be, ga):
    '''
    calculates the volume of a unit cell
    '''
    ca, cb, cg = cos(radians(al)), cos(radians(be)), cos(radians(ga))
    v = a*b*c*sqrt(1+2*ca*cb*cg-ca**2-cb**2-cg**2)
    return v


def dice_coefficient(a, b):
    '''
    dice coefficient 2nt/na + nb
    Compares the similarity of a and b
    :param a: string
    :param b: string
    '''
    a = a.lower()
    b = b.lower()
    if not len(a) or not len(b): return 0.0
    if len(a) == 1:  a=a+'.'
    if len(b) == 1:  b=b+'.'

    a_bigram_list=[]
    for i in range(len(a)-1):
        a_bigram_list.append(a[i:i+2])

    b_bigram_list=[]
    for i in range(len(b)-1):
        b_bigram_list.append(b[i:i+2])

    a_bigrams = set(a_bigram_list)
    b_bigrams = set(b_bigram_list)
    overlap = len(a_bigrams & b_bigrams)
    dice_coeff = overlap * 2.0/(len(a_bigrams) + len(b_bigrams))
    dice_coeff = 1-dice_coeff # invert the result
    if dice_coeff < 0.5:  # make a cutoff for the best matches
        return 0.0
    return round(dice_coeff, 6)


def dice_coefficient2(a,b):
    """
    duplicate bigrams in a word should be counted distinctly
    (per discussion), otherwise 'AA' and 'AAAA' would have a
    dice coefficient of 1...
    """

    if not len(a) or not len(b): return 0.0
    """ quick case for true duplicates """
    if a == b: return 1.0
    """ if a != b, and a or b are single chars, then they can't possibly match """
    if len(a) == 1 or len(b) == 1: return 0.0

    """ use python list comprehension, preferred over list.append() """
    a_bigram_list = [a[i:i+2] for i in range(len(a)-1)]
    b_bigram_list = [b[i:i+2] for i in range(len(b)-1)]

    a_bigram_list.sort()
    b_bigram_list.sort()

    # assignments to save function calls
    lena = len(a_bigram_list)
    lenb = len(b_bigram_list)
    # initialize match counters
    matches = i = j = 0
    while (i < lena and j < lenb):
        if a_bigram_list[i] == b_bigram_list[j]:
            matches += 2
            i += 1
            j += 1
        elif a_bigram_list[i] < b_bigram_list[j]:
            i += 1
        else:
            j += 1

    score = 1-(float(matches)/float(lena + lenb))
    if score < 0.7:
        score = 0.0
    return score


def longest_common_substring(s1, s2):
    m = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]  # @UnusedVariable
    longest, x_longest = 0, 0
    for x in range(1, 1 + len(s1)):
        for y in range(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest - longest: x_longest]


def levenshtein(s1, s2):
    s1 = s1.lower()
    s2 = s2.lower()
    if len(s1) < len(s2):
        return levenshtein(s2, s1)
    if len(s2) == 0:
        return len(s1)
    previous_row = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]

def distance(x1, y1, z1, x2, y2, z2, round_out=False):
    '''
    distance between two points in space
    '''
    import math as m
    d = m.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    if round_out:
        return round(d, round_out)
    else:
        return d

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
    import sys
    v = vol_unitcell(2, 2, 2, 90, 90, 90)
    print(v)
    
    A = ((1, 0, 0),
         (1, 1, 0),
         (0, 0, 1))
    v = (2, 2, 2)
    c = matrix_mult_vector(A, v)
    print(c)
    sys.exit()
    
    from resfile import ResList, ResListEdit
    from atomhandling import FindAtoms
    from dsrparse import DSR_Parser
    import math as m
    
    from dbfile import global_DB
    res_file = 'p21c.res'
    res_list = ResList(res_file)
    reslist = res_list.get_res_list()
    find_atoms = FindAtoms(reslist)
    rle = ResListEdit(reslist, find_atoms)
    dsrp = DSR_Parser(reslist, rle)
    invert = True
    gdb = global_DB(invert)
    db = gdb.build_db_dict()
    fragment = 'PFAnion'
    fragline = gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
    #dbatoms = gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
    dbhead = gdb.get_head_from_fragment(fragment)        # this is only executed once
    resi = True #gdb.get_resi_from_fragment(fragment)
    print(dbhead)
    #uhead = unwrap_head_lines(dbhead)
    #print(uhead)
    # sys.exit()




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
    # sys.exit()
    dsr_string = dsrp.find_dsr_command(line=True).lower()

    regex = 'Q2.*'
    print('found regex in line', find_line(reslist, regex))
    regex2 = '^H[0-9]+\s+'
    multi = find_multi_lines(reslist, regex2)
    print('found multiline in lines', multi)
    if not multi:
        print('nix gefunbden!!!!!!!!')
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

    print('\n\n')
    print(levenshtein('hallo', 'holla'))
    print(dice_coefficient('hallo', 'holla'))
    print(dice_coefficient2('hallo', 'holla'))
    print(which('help'))

    
