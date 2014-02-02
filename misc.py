#-*- encoding: utf-8 -*-
#möp
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

 
def matrix(matrix1,matrix2):
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


def frac_to_cart(frac_coord, cell):
    '''
    Returns an array of catesian coordinates
    '''
    import math as np
    from math import cos, sin, sqrt
    
    cell = [float(y) for y in cell]
    a, b, c, alpha, beta, gamma = cell
    frac_coord = [float(x) for x in frac_coord]
    # convert to radians
    alpha = np.radians(alpha)
    beta  = np.radians(beta)
    gamma = np.radians(gamma)
    
    # cell volume
 
    v = sqrt(1.0 -cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2 + 2*cos(alpha)*cos(beta)*cos(gamma))
 
    tmat = ( [
      [ a,   b*cos(gamma), c*cos(beta)                                       ],
      [ 0.0, b*sin(gamma), c*((cos(alpha)-cos(beta)*cos(gamma))/sin(gamma))  ],
      [ 0.0,          0.0, c*(v/sin(gamma))                                  ]]
      )
 
    cart_coords = matrix([frac_coord], tmat)
    cart_coords = cart_coords[0]
    cart_coords = [str({0:.5}).format(round(i,4)) for i in cart_coords]
    return cart_coords



if __name__ == '__main__':
    from options import OptionsParser 
    from resfile import ResList
    from dsrparse import DSR_Parser
    options = OptionsParser()
    res_list = ResList(options.res_file)
    reslist = res_list.get_res_list()
    dsrp = DSR_Parser(reslist, options)
    
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
    
    coo = frac_to_cart((0.3, 0.3, 0.4), (18, 19, 20, 90, 97, 90))
    #print coo
    print('{:8.5} {:8.5} {:8.5}'.format(*coo), 'neue koordianten')
