# -*- encoding: utf-8 -*-
# 
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

import shutil
from math import cos, sqrt, radians, sin

import os
import random
import re
import string

import mpmath as mpm
from constants import isoatomstr

alphabet = string.ascii_uppercase

__metaclass__ = type  # use new-style classes


def write_file(flist, name):
    """
    Writes the content of list to name.
    :param flist: list
    :param name:  string
    """
    with open(name, 'w') as ofile:
        for line in flist:  # modified reslist
            ofile.write("%s" % line)  # write the new file


def extract_tarfile(file, targetdir):
    """
    Extracts .tar.gz "file" to "targetdir"
    :param file: .tar.gz file
    :param targetdir: target directory
    :return: bool for sucess
    """
    import tarfile
    try:
        with tarfile.open(file) as tarobj:
            tarobj.extractall(path=targetdir)
    except tarfile.ReadError as e:
        print('*** Could not extract tarfile ***')
        print('***', e, '***')
        return False
    return True


def join_floats(float_list, places=3, ):
    """
    >>> l = [1.23456789123456, 1, 2, 3, 3.45, 5.6543, 1,3456]
    >>> join_floats(l)
    '1.235 1.000 2.000 3.000 3.450 5.654 1.000 3456.000'
    """
    str_floats = " ".join(format(i, "{}.{}f".format(places + 2, places)) for i in float_list)
    return str_floats


def touch(fname, times=None):
    """
    creates an empty file like unix touch command
    """
    with open(fname, 'a'):
        os.utime(fname, times)


def check_file_exist(filename):
    """
    Check if a file exists and has some content.
    A file zize above 0 returns True.
    A non-existing file returns False.
    A file size of 0 retrurns 'zero'.

    # keep these paths for the unit tests!
    >>> check_file_exist('p21c.res')
    True
    >>> check_file_exist('foo.bar')
    File "foo.bar" not found!
    False
    >>> check_file_exist('empty.txt')
    'zero'
    >>> check_file_exist('../misc.py')
    True
    """
    status = False
    if os.path.isfile(filename):
        filesize = int(os.stat(str(filename)).st_size)
    else:
        print('File "{}" not found!'.format(filename))
        return False
    if isinstance(filesize, int) and filesize > 0:
        status = True
    if isinstance(filesize, int) and filesize == 0:
        status = 'zero'
    return status


def walkdir(rootdir, include="", exclude=None):
    """
    Returns a list of files in all subdirectories with full path.
    :param rootdir: base path from which walk should start
    :param include: list of file endings to include only e.g. ['.py', '.res']
    :return: list of files

    >>> walkdir("../docs") #doctest: +REPORT_NDIFF +NORMALIZE_WHITESPACE +ELLIPSIS
    ['../docs/test.txt']
    >>> walkdir("../setup/modpath.iss")
    ['../setup/modpath.iss']
    >>> walkdir("../setup/modpath.iss", exclude=['.iss'])
    []
    >>> walkdir("../docs", exclude=['.txt']) #doctest: +REPORT_NDIFF +NORMALIZE_WHITESPACE +ELLIPSIS
    []
    """
    if exclude is None:
        exclude = [""]
    results = []
    if not os.path.isdir(rootdir):
        if os.path.splitext(rootdir)[1] in exclude:
            return []
        return [rootdir]
    for root, subFolders, files in os.walk(rootdir):
        for file in files:
            fullfilepath = os.path.join(root, file)
            if exclude:
                if os.path.splitext(fullfilepath)[1] in exclude:
                    continue
            if include:
                if os.path.splitext(fullfilepath)[1] in include:
                    results.append(os.path.normpath(fullfilepath).replace('\\', '/'))
            else:
                results.append(os.path.normpath(fullfilepath).replace('\\', '/'))
    return results


def pairwise(iterable):
    """
     s -> (s0,s1), (s2,s3), (s4, s5), ...

     >>> liste = ['C1', 'C2', 'C2', 'C3', 'C4', 'C5', 'C5', 'C6']
     >>> pairwise(liste)
     [('C1', 'C2'), ('C2', 'C3'), ('C4', 'C5'), ('C5', 'C6')]
     """
    a = iter(iterable)
    return list(zip(a, a))


def mean(values):
    """
    returns mean value of a list of numbers

    >>> mean([1, 2, 3, 4, 1, 2, 3, 4])
    2.5
    >>> round(mean([1, 2, 3, 4, 1, 2, 3, 4.1, 1000000]), 4)
    111113.3444
    """
    mean = sum(values) / float(len(values))
    return mean


def median(nums):
    """
    calculates the median of a list of numbers
    >>> median([2])
    2
    >>> median([1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4])
    2.5
    >>> median([1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4.1, 1000000])
    3
    >>> median([])
    Traceback (most recent call last):
    ...
    ValueError: Need a non-empty iterable
    """
    ls = sorted(nums)
    n = len(ls)
    if n == 0:
        raise ValueError("Need a non-empty iterable")
    # for uneven list length:
    elif n % 2 == 1:
        # // is floordiv:
        return ls[n // 2]
    else:
        return sum(ls[int(int(n) / 2 - 1):int(int(n) / 2 + 1)]) / 2.0


def std_dev(data):
    """
    returns standard deviation of values rounded to pl decimal places
    S = sqrt( (sum(x-xm)^2) / n-1 )
    xm = sum(x)/n
    >>> l1 = [1.334, 1.322, 1.345, 1.451, 1.000, 1.434, 1.321, 1.322]
    >>> l2 = [1.234, 1.222, 1.345, 1.451, 2.500, 1.234, 1.321, 1.222]
    >>> round(std_dev(l1), 8)
    0.13797871
    >>> round(std_dev(l2), 8)
    0.43536797
    >>> median(l1)
    1.328
    >>> mean(l1)
    1.316125
    """
    if len(data) == 0:
        return 0
    K = data[0]
    n = 0
    Sum = 0
    Sum_sqr = 0
    for x in data:
        n += 1
        Sum += x - K
        Sum_sqr += (x - K) * (x - K)
    variance = (Sum_sqr - (Sum * Sum) / n) / (n - 1)
    # use n instead of (n-1) if want to compute the exact variance of the given data
    # use (n-1) if data are samples of a larger population
    return sqrt(variance)


def nalimov_test(data):
    """
    returns a index list of outliers base on the Nalimov test for data.
    Modified implementation of:
    "R. Kaiser, G. Gottschalk, Elementare Tests zur Beurteilung von Messdaten
    Bibliographisches Institut, Mannheim 1972." 
    
    >>> data = [1.120, 1.234, 1.224, 1.469, 1.145, 1.222, 1.123, 1.223, 1.2654, 1.221, 1.215]
    >>> nalimov_test(data)
    [3]
    """
    # q-values for degrees of freedom:
    f = {1 : 1.409, 2: 1.645, 3: 1.757, 4: 1.814, 5: 1.848, 6: 1.870, 7: 1.885, 8: 1.895,
         9 : 1.903, 10: 1.910, 11: 1.916, 12: 1.920, 13: 1.923, 14: 1.926, 15: 1.928,
         16: 1.931, 17: 1.933, 18: 1.935, 19: 1.936, 20: 1.937, 30: 1.945}
    fact = sqrt(float(len(data)) / (len(data) - 1))
    fval = len(data) - 2
    if fval < 2:
        return []
    outliers = []
    if fval in f:
        # less strict than the original:
        q_crit = f[fval]
    else:
        q_crit = 1.95
    for num, i in enumerate(data):
        q = abs(((i - median(data)) / std_dev(data)) * fact)
        if q > q_crit:
            outliers.append(num)
    return outliers


def flatten(lis):
    """
    Given a list, possibly nested to any level, return it flattened.
    From: http://code.activestate.com/recipes/578948-flattening-an-arbitrarily-nested-list-in-python/

    >>> flatten([['wer', 234, 'brdt5'], ['dfg'], [[21, 34,5], ['fhg', 4]]])
    ['wer', 234, 'brdt5', 'dfg', 21, 34, 5, 'fhg', 4]
    """
    new_lis = []
    for item in lis:
        if type(item) == type([]):
            new_lis.extend(flatten(item))
        else:
            new_lis.append(item)
    return new_lis


def sortedlistdir(directory):
    """
    returns a sorted list of files in directory directory.
    :param directory: directory
    :type directory: string
    >>> sortedlistdir("../old")
    ['dsr.py']
    >>> sortedlistdir("foobar/")
    False
    """
    try:
        dirlist = os.listdir(directory)
    except(IOError, OSError):
        return False
    dirlist.sort()
    return dirlist


def find_line_of_residue(reslist, resinumber):
    """
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
    >>> find_line_of_residue(['REM TEST', 'C1 1 -0.00146 0.26814 0.06351 11.00 0.05', \
        'RESI 4 BENZ', 'C2 1 -1.13341 -0.23247 -0.90730 11.00 0.05'], "4")
    [2, 'RESI 4 BENZ']
    """
    for n, line in enumerate(reslist):
        if line.upper().startswith('RESI'):
            if line.split()[1] == str(resinumber):
                return [n, line]


def ll_to_string(inputlist):
    """
    converts list of list to string with four whitespaces between each list element
    >>> ll_to_string([['REM', 'TEST', 'C1'], ['1', '-0.00146', '0.26814', '0.06351']])
    'REM   TEST   C1\\n1   -0.00146   0.26814   0.06351'
    """
    inputlist = [list(map(str, i)) for i in inputlist]
    newlist = []
    for i in inputlist:
        newlist.append('   '.join(i).rstrip())
    string = '\n'.join(newlist)
    return string


def multiline_test(line):
    """
    test if the current line is a multiline with "=" at the end
    :param line: 'O1 3 -0.01453 1.66590 0.10966 11.00 0.05 ='
    :type line: string
    >>> line = 'C1    1    0.278062    0.552051    0.832431    11.00000    0.02895    0.02285 ='
    >>> multiline_test(line)
    True
    >>> line = 'C1    1    0.278062    0.552051    0.832431    11.00000    0.05 '
    >>> multiline_test(line)
    False
    """
    line = line.rpartition('=')  # partition the line in before, sep, after
    line = ''.join(line[0:2])  # use all including separator
    line = line.rstrip()  # strip spaces
    if line.endswith('='):
        return True
    else:
        return False


def find_line(inputlist, regex, start=None):
    """
    returns the index number of the line where regex is found in the inputlist
    if stop is true, stop searching with first line found
    :param inputlist: list of strings
    :type inputlist: list
    :param regex: regular expression to search
    :type regex: string
    :param start: line number where to start the search
    :param start: start searching at line start
    :type start: string or int
    >>> inp = ['Hallo blub', 'foo bar blub', '123', '1 blub 2 3 4']
    >>> find_line(inp, '.*blub.*')
    0
    >>> find_line(inp, 'nonono')
    False
    >>> inp = [['foo'],['bar']]
    >>> find_line(inp, '.*blub.*') #doctest: +REPORT_NDIFF +NORMALIZE_WHITESPACE +ELLIPSIS
    Traceback (most recent call last):
        ...
    TypeError: expected string or ...
    """
    if start:
        start = int(start)
        inputlist_slice = inputlist[start:]
        for i, string in enumerate(inputlist_slice, start):
            if re.match(regex, string, re.IGNORECASE):
                return i  # returns the index number if regex found
    else:
        for i, string in enumerate(inputlist):
            if re.match(regex, string, re.IGNORECASE):
                return i  # returns the index number if regex found
    return False  # returns False if no regex found (xt solution has no fvar)


def remove_line(reslist, linenum, rem=False, remove=False, frontspace=False):
    """
    removes a single line from the res file with tree different methods.
    The default is a space character in front of the line (frontspace).
    This removes the line in the next refinement cycle. "rem" writes rem
    in front of the line and "remove" clears the line.
    :param reslist: .res file list
    :param linenum: integer, line number
    :param rem:     True/False, activate comment with 'REM' in front
    :param remove:  True/False, remove the line
    :param frontspace: True/False, activate removing with a front space
    """
    line = reslist[linenum]
    if rem:  # comment out with 'rem ' in front
        reslist[linenum] = 'rem ' + line
        if multiline_test(line):
            reslist[linenum + 1] = 'rem ' + reslist[linenum + 1]
    elif remove:  # really delete the line "linenum"
        if multiline_test(line):
            reslist[linenum] = ''
            reslist[linenum + 1] = ''
        else:
            reslist[linenum] = ''
    if frontspace:  # only put a space in front
        reslist[linenum] = ' ' + line
        if multiline_test(line):
            reslist[linenum + 1] = ' ' + reslist[linenum + 1]
    return reslist


def find_multi_lines(inputlist, regex):
    """
    returns the index number of all lines where regex is found in the inputlist
    ! this method is case insensitive !
    >>> input = ['Hallo blub', 'foo bar blub', '123', '1 blub 2 3 4']
    >>> find_multi_lines(input, '.*blub.*')
    [0, 1, 3]
    >>> find_multi_lines(input, 'blabla')
    []
    >>> inp = [['foo'],['bar']]
    >>> find_multi_lines(inp, r'.*blub.*') # doctest: +NORMALIZE_WHITESPACE +REPORT_NDIFF +ELLIPSIS
    Traceback (most recent call last):
        ...
    TypeError: expected string or ...
    """
    reg = re.compile(regex, re.IGNORECASE)
    foundlist = []
    for i, string in enumerate(inputlist):
        if reg.match(string):
            foundlist.append(i)  # returns the index number if regex found
    return foundlist


def remove_file(filename, exit_dsr=False, terminate=False):
    """
    removes the file "filename" from disk
    program exits when exit is true
    platon gets terminated if terminate is true

    >>> remove_file('foobar')
    """
    if os.path.isfile(filename):
        try:
            os.remove(filename)
        except(IOError, OSError):
            print('can not delete {}'.format(filename))
            # print 'unable to cleanup ins {} files!'.format(file)
            if terminate:
                pgrogname = terminate
                pgrogname.terminate()
            if exit_dsr:
                sys.exit()


def copy_file(source, target):
    """
    Copy a file from source to target. Source can be a single file or
    a directory. Target can be a single file or a directory.
    :param source: list or string
    :param target: string
    """
    target_path = os.path.dirname(target)
    source_file = os.path.basename(source)
    listcopy = False
    if isinstance(source, (list, tuple)):
        listcopy = True
    if not os.path.exists(target_path) and target_path != '':
        try:
            os.makedirs(target_path)
        except(IOError, OSError):
            print('Unable to create directory {}.'.format(target_path))
    try:
        if listcopy:
            for filen in source:
                shutil.copy(filen, target)
        else:
            shutil.copy(source, target)
    except IOError as e:
        print('Unable to copy {}.'.format(source_file))
        print(e)


def make_directory(dirpath):
    """
    create a directory with all subdirs from the last existing path
    :param dirpath: string
    """
    try:
        os.makedirs(dirpath)
    except(IOError, OSError):
        print('Unable to create directory {}.'.format(dirpath))
        return False
    return True


def wrap_headlines(dbhead, width=75):
    """
    wraps lines of a restraint header to prevent too long lines in
    SHELXL. wrapping is done with = at the end of a line and ' ' at
    start of the next line
    :param dbhead: header with restraints
    :param width: wrap after width characters
    >>> line = ['foo bar this is text to wrap. blah bub']
    >>> wrap_headlines(line, 15)
    ['foo bar this is =\\n   text to wrap. =\\n   blah bub\\n']
    >>> wrap_headlines(['SADI C1 C2 C3 C4'], 10)
    ['SADI C1 C2 =\\n   C3 C4\\n']
    """
    import textwrap
    for num, line in enumerate(dbhead):
        line = textwrap.wrap(line, width, subsequent_indent='  ')
        if len(line) > 1:
            newline = []
            for n, l in enumerate(line):
                if n < len(line) - 1:
                    l += ' =\n'
                newline.append(l)
            dbhead[num] = ' '.join(newline)
    for num, line in enumerate(dbhead):
        line = ' '.join(line.strip().split(' '))
        dbhead[num] = line + '\n'
    return dbhead


def wrap_stringlist(strlist, width=75):
    """
    Wraps the text lines of a list to width characters.
    >>> wrap_stringlist(['REM foo bar baz foo bar baz blah blub', 'REM 2foo bar baz foo bar baz blah blub'], 36)
    ['REM foo bar baz foo bar baz blah\\nREM blub\\n', 'REM 2foo bar baz foo bar baz blah\\nREM blub\\n']
    """
    import textwrap
    wrapped = []
    for num, line in enumerate(strlist):
        wrapped.append('\n'.join(textwrap.wrap(line, width, subsequent_indent='REM ')) + '\n')
    return wrapped


def wrap_text(inText, maxlen=70, subsequent_indent='=\n'):
    """
    Text wrapper without need for textwrap package.
    >>> wrap_text('SADI 0.02 C1A C2A C2A C3A C3A C4A C4A C5A C5A C6A', maxlen=30)
    'SADI 0.02 C1A C2A C2A C3A C3A =\\n C4A C4A C5A C5A C6A'
    """
    line_list = []
    wrapped = []
    inText_list = inText.strip().split()
    for n, word in enumerate(inText_list):
        testline = ' '.join(line_list)
        if len(word + testline + subsequent_indent) >= maxlen:
            if n < len(inText_list) - 1:
                line_list.append(word + ' ' + subsequent_indent)
            else:
                line_list.append(word)
            wrapped.append(' '.join(line_list).strip(' '))
            line_list = []
        else:
            line_list.append(word)
    if line_list:
        wrapped.append(' '.join(line_list).strip(' '))
    return ' '.join(wrapped)


def unwrap_head_lines(headlines):
    """
    if a line is wrapped like "SADI C1 C2 =\n", "  C3 C4" or "SADI C1 C2=\n", "  C3 C4"
    this function returns "SADI C1 C2 C3 C4"
    :type headlines: list
    :param headlines: list of strings from a res file
    >>> unwrap_head_lines(["SADI C1 C2 =\\n", "  C3 C4"])
    ['SADI C1 C2 C3 C4']
    >>> unwrap_head_lines(['foo bar this is =\\n   text to wrap. =\\n   blah bub\\n'])
    ['foo bar this is text to wrap. blah bub']
    >>> unwrap_head_lines(wrap_headlines(["REM name: ser\\nREM HFIX 23 C2 C3 C4  C6 C8 C10\\nSADI bar"]))
    ['REM name: ser\\nREM HFIX 23 C2 C3 C4  C6 C8 C10\\nSADI bar\\n']
    """
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
        line = line.strip(' \r\n=').replace('=', ' ')
        tmp = tmp + ' ' + line
    spline = tmp.split()
    last = ''
    for n, i in enumerate(spline):
        command = i[:4].upper()
        # leave out 'REM COMMAND':
        if command in constants.SHX_CARDS and last != 'REM':
            spline[n] = '\n' + spline[n]
        last = command
    new_head = ' '.join(spline).strip().split('\n')
    new_head = [i.rstrip(' ') for i in new_head]
    return new_head


def which(name, flags=os.X_OK, exts=None):
    """
    Search PATH for executable files with the given name.

    On MS-Windows the only flag that has any meaning is os.F_OK. Any other
    flags will be ignored.
    #>>> which('dsr.bat') # doctest: +NORMALIZE_WHITESPACE +REPORT_NDIFF +ELLIPSIS
    ['d:\\\Programme\\\DSR\\\dsr.bat']
    """
    if exts is None:
        exts = ['.exe', '.EXE', '.bat']
    result = []
    # exts = filter(None, os.environ.get('PATHEXT', '').split(os.pathsep))
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
    """
    strips the part symbol like C1_4b from an atom name
    :param atom: 'C1_4b'
    :type atom: string

    >>> remove_partsymbol('C2_4b')
    'C2_4'
    >>> remove_partsymbol('C22_b')
    'C22'
    >>> remove_partsymbol('C_5')
    'C_5'
    >>> remove_partsymbol('SAME/SADI')
    'SAME/SADI'
    >>> remove_partsymbol('C22_4^b')
    'C22_4'
    >>> remove_partsymbol('C23^b')
    'C23'
    >>> remove_partsymbol('C24_0^b')
    'C24'
    >>> remove_partsymbol('C25_0b')
    'C25'
    """
    if '_' in atom or '^' in atom:
        if '^' in atom:
            # since SHELXL 2016/5, residue and part are devided by "^"
            name = atom.split('^')[0]
            if name.split('_')[-1] == '0':
                return name.split('_')[0]
            else:
                return name
        else:
            presuff = atom.split('_')
            prefix, suffix = presuff[0], presuff[-1].strip(string.ascii_letters)
        if not suffix:
            return prefix
        else:
            if suffix == '0':
                atom = prefix
            else:
                atom = prefix + '_' + suffix
    return atom


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    """
    returns a random ID like 'L5J74W'
    :param size: length of the string
    :type size: integer
    :param chars: characters used for the ID
    :type chars: string

    >>> id_generator(1, 'a')
    'a'
    """
    return ''.join(random.choice(chars) for _ in range(size))


def shift(seq, n):
    """
    left-shift a sliceable object by n
    :param seq: sequence to shift
    :type seq: string or list
    :param n: shift length
    :type n: int
    >>> shift('hello world', 3)
    'lo worldhel'
    >>> shift(['sdfg', 'dsfg', '111', '222'], 1)
    ['dsfg', '111', '222', 'sdfg']
    """
    n %= len(seq)
    return seq[n:] + seq[:n]


def atomic_distance(p1, p2, cell, shortest_dist=False):
    """
    p1 and p2 are x, y , z coordinates as list ['x', 'y', 'z']
    cell are the cell parameters as list: ['a', 'b', 'c', 'alpha', 'beta', 'gamma']

    Returns the distance between the two points (Atoms). If shortest_dist is True, the
    shortest distance ignoring translation is computed.

    >>> cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
    >>> coord1 = (-0.186843,   0.282708,   0.526803)
    >>> coord2 = (-0.155278,   0.264593,   0.600644)
    >>> atomic_distance(coord1, coord2, cell)
    1.5729229943265979
    >>> cell = [14.0452, 14.0452, 12.2008, 90.000, 90.000, 90.000]
    >>> atomic_distance([-0.184235, 0.955143, -0.181577], [0.437437, 1.162938, 0.182932], cell)
    10.224259725884533
    >>> atomic_distance([0.455143, 1.184235, 0.181577], [0.437437, 1.162938, 0.182932], cell, shortest_dist=True)
    0.38934604708396064
    """
    a, b, c = cell[:3]
    al = radians(cell[3])
    be = radians(cell[4])
    ga = radians(cell[5])
    if shortest_dist:
        x1, y1, z1 = [x + 99.5 for x in p1]
        x2, y2, z2 = p2
        dx = (x1 - x2) % 1 - 0.5
        dy = (y1 - y2) % 1 - 0.5
        dz = (z1 - z2) % 1 - 0.5
    else:
        x1, y1, z1 = p1
        x2, y2, z2 = p2
        dx = (x1 - x2)
        dy = (y1 - y2)
        dz = (z1 - z2)
    dsq = (a * dx) ** 2 + (b * dy) ** 2 + (c * dz) ** 2 + 2 * b * c * cos(al) * dy * dz + \
          2 * dx * dz * a * c * cos(be) + 2 * dx * dy * a * b * cos(ga)
    return sqrt(dsq)


def frac_to_cart(frac_coord, cell):
    """
    Converts fractional coordinates to cartesian coodinates
    :param frac_coord: [float, float, float]
    :param cell:       [float, float, float, float, float, float]
    >>> import mpmath as mpm
    >>> cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
    >>> coord1 = (-0.186843,   0.282708,   0.526803)
    >>> print(frac_to_cart(coord1, cell))
    [-2.741505423999065, 5.909586678000002, 10.775200700893734]
    """
    a, b, c, alpha, beta, gamma = cell
    x, y, z = frac_coord
    alpha = radians(alpha)
    beta = radians(beta)
    gamma = radians(gamma)
    try:
        cosastar = (cos(beta) * cos(gamma) - cos(alpha)) / (sin(beta) * sin(gamma))
        sinastar = sqrt(1 - cosastar ** 2)
    except ValueError:
        print("*** Malformed unit cell parameters found! Please correct the database entry. ***")
        sys.exit()
    Xc = a * x + (b * cos(gamma)) * y + (c * cos(beta)) * z
    Yc = 0 + (b * sin(gamma)) * y + (-c * sin(beta) * cosastar) * z
    Zc = 0 + 0 + (c * sin(beta) * sinastar) * z
    return [Xc, Yc, Zc]


class A(object):
    """
    orthogonalization matrix
    e.g. converts fractional coordinates to cartesian coodinates
    >>> import mpmath as mpm
    >>> cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
    >>> coord = (-0.186843,   0.282708,   0.526803)
    >>> A = A(cell).orthogonal_matrix
    >>> print(mpm.nstr(A*mpm.matrix(coord)))
    [-2.74151]
    [ 5.90959]
    [ 10.7752]
    >>> cartcoord = mpm.matrix([['-2.74150542399906'], ['5.909586678'], ['10.7752007008937']])
    >>> print(mpm.nstr(A**-1*cartcoord))
    [-0.186843]
    [ 0.282708]
    [ 0.526803]
    """

    def __init__(self, cell):
        self.a, self.b, self.c, alpha, beta, gamma = cell
        self.V = vol_unitcell(self.a, self.b, self.c, alpha, beta, gamma)
        self.alpha = radians(alpha)
        self.beta = radians(beta)
        self.gamma = radians(gamma)

    @property
    def orthogonal_matrix(self):
        """
        Converts von fractional to cartesian.
        Invert the matrix to do the opposite.
        """
        import mpmath as mpm
        Am = mpm.matrix([[self.a, self.b * cos(self.gamma), self.c * cos(self.beta)],
                         [0, self.b * sin(self.gamma),
                          (self.c * (cos(self.alpha) - cos(self.beta) * cos(self.gamma)) / sin(self.gamma))],
                         [0, 0, self.V / (self.a * self.b * sin(self.gamma))]])
        return Am


def cart_to_frac(cart_coord, cell):
    """
    converts cartesian coordinates to fractional coordinates
    :param cart_coord: [float, float, float]
    :param cell:       [float, float, float, float, float, float]
    >>> cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
    >>> coords = [-2.74150542399906, 5.909586678, 10.7752007008937]
    >>> cart_to_frac(coords, cell)
    [-0.186843, 0.282708, 0.526803]
    """
    a, b, c, alpha, beta, gamma = cell
    X, Y, Z = cart_coord
    alpha = radians(alpha)
    beta = radians(beta)
    gamma = radians(gamma)
    cosastar = (cos(beta) * cos(gamma) - cos(alpha)) / (sin(beta) * sin(gamma))
    sinastar = sqrt(1 - cosastar ** 2)
    z = Z / (c * sin(beta) * sinastar)
    y = (Y - (-c * sin(beta) * cosastar) * z) / (b * sin(gamma))
    x = (X - (b * cos(gamma)) * y - (c * cos(beta)) * z) / a
    return [round(x, 8), round(y, 8), round(z, 8)]


def zero(m, n):
    """
    Create zero matrix of dimension m,n
    :param m: integer
    :param n: integer

    >>> zero(5, 3)
    [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]
    """
    new_matrix = [[0 for row in range(n)] for col in range(m)]  # @UnusedVariable
    return new_matrix


def determinante(a):
    """
    return determinant of 3x3 matrix
    Deprecated, use mpmath instead!!!

    >>> m1 = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    >>> determinante(m1)
    8
    """
    return (a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2])
            - a[1][0] * (a[0][1] * a[2][2] - a[2][1] * a[0][2])
            + a[2][0] * (a[0][1] * a[1][2] - a[1][1] * a[0][2]))


def subtract_vect(a, b):
    """
    subtract vector b from vector a
    Deprecated, use mpmath instead!!!
    :param a: [float, float, float]
    :param b: [float, float, float]

    >>> subtract_vect([1, 2, 3], [3, 2, 2])
    (-2, 0, 1)
    """
    return (a[0] - b[0],
            a[1] - b[1],
            a[2] - b[2])


def matrix_minus_vect(m, v):
    """
    >>> source = [[-0.01453, 1.6659, 0.10966], [-0.00146, 0.26814, 0.06351], [-0.27813, -0.21605, 1.52795]]
    >>> cent = [-0.09804,     0.57266333 , 0.56704]
    >>> matrix_minus_vect(source, cent) # DOCTEST: +NORMALIZE_WHITESPACE +ELLIPSIS
    [[0.08351, 1.09323667, -0.45738], [0.09658, -0.30452333000000004, -0.50353], [-0.18008999999999997, -0.78871333, 0.9609099999999999]]
    >>> matrix_minus_vect([[1,1,1],[2,2,2],[3,3,3]], [1,1,1])
    [[0, 0, 0], [1, 1, 1], [2, 2, 2]]
    """
    result = []
    for coords in m:
        result.append([coords[0] - v[0], coords[1] - v[1], coords[2] - v[2]])
    return result


def matrix_plus_vect(m, v):
    """
    >>> source = [[-0.01453, 1.6659, 0.10966], [-0.00146, 0.26814, 0.06351], [-0.27813, -0.21605, 1.52795]]
    >>> cent = [-0.09804,     0.57266333 , 0.56704]
    >>> matrix_plus_vect(source, cent)
    [[-0.11257, 2.23856333, 0.6767], [-0.0995, 0.84080333, 0.6305499999999999], [-0.37617, 0.35661333000000006, 2.09499]]
    >>> matrix_minus_vect([[1,1,1],[2,2,2],[3,3,3]], [1,1,1])
    [[0, 0, 0], [1, 1, 1], [2, 2, 2]]
    """
    result = []
    for coords in m:
        result.append([coords[0] + v[0], coords[1] + v[1], coords[2] + v[2]])
    return result


def transpose(a):
    """
    transposes a matrix

    >>> m = [[1, 2, 3], [1, 2, 3], [1, 2, 3]]
    >>> transpose(m)
    [(1, 1, 1), (2, 2, 2), (3, 3, 3)]
    """
    return list(zip(*a))


def norm_vec(a):
    """
    returns a normalized vector

    >>> norm_vec([1, 2, 1])
    (0.4082482904638631, 0.8164965809277261, 0.4082482904638631)
    """
    l = sqrt(a[0] ** 2 + a[1] ** 2 + a[2] ** 2)
    return a[0] / l, a[1] / l, a[2] / l


def vol_tetrahedron(a, b, c, d, cell=None):
    """
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

    >>> cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
    >>> a = (0.838817,   0.484526,   0.190081) # a ist um 0.01 ausgelenkt
    >>> b = (0.875251,   0.478410,   0.256955)
    >>> c = (0.789290,   0.456520,   0.301616)
    >>> d = (0.674054,   0.430194,   0.280727)
    >>> print('volume of Benzene ring atoms:')
    volume of Benzene ring atoms:
    >>> print(round(vol_tetrahedron(a, b, c, d, cell), 7))
    0.0633528
    """
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
    volume = abs((D / 6))
    return volume


def vol_unitcell(a, b, c, al, be, ga):
    """
    calculates the volume of a unit cell
    >>> v = vol_unitcell(2, 2, 2, 90, 90, 90)
    >>> print(v)
    8.0
    """
    ca, cb, cg = cos(radians(al)), cos(radians(be)), cos(radians(ga))
    v = a * b * c * sqrt(1 + 2 * ca * cb * cg - ca ** 2 - cb ** 2 - cg ** 2)
    return v


def dice_coefficient(a, b, case_insens=True):
    """
    :type a: str
    :type b: str
    :type case_insens: bool
    dice coefficient 2nt/na + nb.
    https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Dice%27s_coefficient#Python
    >>> dice_coefficient('hallo', 'holla')
    0.25
    >>> dice_coefficient('Banze', 'Benzene')
    0.444444
    >>> dice_coefficient('halo', 'Haaallo')
    0.75
    >>> dice_coefficient('hallo', 'Haaallo')
    0.888889
    >>> dice_coefficient('hallo', 'Hallo')
    1.0
    >>> dice_coefficient('aaa', 'BBBBB')
    0.0
    """
    if case_insens:
        a = a.lower()
        b = b.lower()
    if not len(a) or not len(b):
        return 0.0
    if len(a) == 1:
        a = a + u'.'
    if len(b) == 1:
        b = b + u'.'
    a_bigram_list = []
    for i in range(len(a) - 1):
        a_bigram_list.append(a[i:i + 2])
    b_bigram_list = []
    for i in range(len(b) - 1):
        b_bigram_list.append(b[i:i + 2])
    a_bigrams = set(a_bigram_list)
    b_bigrams = set(b_bigram_list)
    overlap = len(a_bigrams & b_bigrams)
    dice_coeff = overlap * 2.0 / (len(a_bigrams) + len(b_bigrams))
    return round(dice_coeff, 6)


def dice_coefficient2(a, b, case_insens=True):
    """ 
    :type a: str
    :type b: str
    :type case_insens: bool
    duplicate bigrams in a word should be counted distinctly
    (per discussion), otherwise 'AA' and 'AAAA' would have a 
    dice coefficient of 1...
    https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Dice%27s_coefficient#Python

    This implementation is reverse. 1 means not hit, 0 means best match
    >>> dice_coefficient2('hallo', 'holla')
    0.5
    >>> dice_coefficient2('Banze', 'Benzene')
    0.8
    >>> dice_coefficient2('halo', 'Haaallo')
    1.333333
    >>> dice_coefficient2('hallo', 'Haaallo')
    1.6
    >>> dice_coefficient2('hallo', 'Hallo')
    2.0
    >>> dice_coefficient2('aaa', 'BBBBB')
    0.0
    >>> dice_coefficient2('', '')
    0.0
    """
    if case_insens:
        a = a.lower()
        b = b.lower()
    if not len(a) or not len(b):
        return 0.0
    # quick case for true duplicates
    if a == b:
        return 2.0
    # if a != b, and a or b are single chars, then they can't possibly match
    if len(a) == 1 or len(b) == 1:
        return 0.0
    # use python list comprehension, preferred over list.append()
    a_bigram_list = [a[i:i + 2] for i in range(len(a) - 1)]
    b_bigram_list = [b[i:i + 2] for i in range(len(b) - 1)]
    a_bigram_list.sort()
    b_bigram_list.sort()
    # assignments to save function calls
    lena = len(a_bigram_list)
    lenb = len(b_bigram_list)
    # initialize match counters
    matches = i = j = 0
    while i < lena and j < lenb:
        if a_bigram_list[i] == b_bigram_list[j]:
            matches += 2
            i += 1
            j += 1
        elif a_bigram_list[i] < b_bigram_list[j]:
            i += 1
        else:
            j += 1
    score = float(2*matches) / float(lena + lenb)
    return round(score, 6)


def longest_common_substring(s1, s2):
    """
    returns the longest common substring of two strings
    :param s1: a string
    :type s1: str
    :param s2: a second string
    :type s2: str

    >>> longest_common_substring('hello world how is foo bar?', 'hello daniel how is foo in the world?')
    ' how is foo '
    """
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
    """
    >>> levenshtein('hallo', 'holla')
    2
    >>> dice_coefficient('hallo', 'holla')
    0.25
    """
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
            # j+1 instead of j since previous_row and current_row are one character longer:
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1  # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]


def distance(x1, y1, z1, x2, y2, z2, round_out=0):
    """
    distance between two points in space for orthogonal axes.
    >>> distance(1, 1, 1, 2, 2, 2, 4)
    1.7321
    >>> distance(1, 0, 0, 2, 0, 0, 4)
    1.0
    """
    from math import sqrt
    d = sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
    if round_out:
        return round(d, round_out)
    else:
        return d


def coord_to_shx_atom(coordinates):
    """
    Transformes a list of coordinates to a list of shelxl atom strings.
    """
    strlist = []
    sfac_num = 1
    for num, coord in enumerate(coordinates):
        at = "AX" + str(num)
        s = isoatomstr.format(at, sfac_num, coord[0], coord[1], coord[2], 11.0000, 0.03)
        strlist.append(s)
    return strlist


def calc_ellipsoid_axes(coords, uvals, cell, probability=0.5, longest=True):
    """
    This method calculates the principal axes of an ellipsoid as list of two
    fractional coordinate triples.
    Many thanks to R. W. Grosse-Kunstleve and P. D. Adams
    for their great publication on the handling of atomic anisotropic displacement
    parameters:
    R. W. Grosse-Kunstleve, P. D. Adams, J Appl Crystallogr 2002, 35, 477–480.

    F = ... * exp ( -2π²[ h²(a*)²U11 + k²(b*)²U22 + ... + 2hka*b*U12 ] )

    SHELXL atom:
    Name type  x      y      z    occ     U11 U22 U33 U23 U13 U12
    F3    4    0.210835   0.104067   0.437922  21.00000   0.07243   0.03058 =
       0.03216  -0.01057  -0.01708   0.03014
    >>> import mpmath as mpm
    >>> cell = [10.5086, 20.9035, 20.5072, 90, 94.13, 90]
    >>> coords = [0.210835,   0.104067,   0.437922]
    >>> uvals = [0.07243, 0.03058, 0.03216, -0.01057, -0.01708, 0.03014]
    >>> l = calc_ellipsoid_axes(coords, uvals, cell, longest=True)
    >>> print(mpm.nstr(l))
    [[0.24765096, 0.11383281, 0.43064756], [0.17401904, 0.09430119, 0.44519644]]
    >>> calc_ellipsoid_axes(coords, uvals, cell, longest=False)
    [[[0.24765096, 0.11383281, 0.43064756], [0.218406, 0.09626142, 0.43746127], [0.21924358, 0.10514684, 0.44886868]], [[0.17401904, 0.09430119, 0.44519644], [0.203264, 0.11187258, 0.43838273], [0.20242642, 0.10298716, 0.42697532]]]
    >>> cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
    >>> coords = [0.210835,   0.104067,   0.437922]
    >>> uvals = [0.07243, -0.03058, 0.03216, -0.01057, -0.01708, 0.03014]
    >>> calc_ellipsoid_axes(coords, uvals, cell, longest=True)
    <BLANKLINE>
    Ellipsoid is non positive definite!
    <BLANKLINE>
    False

    >>> uvals = [0.07243, 0.03058, 0.03216, -0.01057, -0.01708]
    >>> calc_ellipsoid_axes(coords, uvals, cell, longest=False)
    Traceback (most recent call last):
    ...
    Exception: 6 Uij values have to be supplied!

    >>> cell = (10.5086, 20.9035, 90, 94.13, 90)
    >>> coords = [0.210835,   0.104067,   0.437922]
    >>> uvals = [0.07243, 0.03058, 0.03216, -0.01057, -0.01708, 0.03014]
    >>> calc_ellipsoid_axes(coords, uvals, cell, longest=True)
    Traceback (most recent call last):
    ...
    Exception: cell needs six parameters!

    :param coords: coordinates of the respective atom in fractional coordinates
    :type coords: list
    :param uvals: Uij valiues of the respective ellipsoid on fractional
                  basis like in cif and SHELXL format
    :type uvals: list
    :param cell: unit cell of the structure: a, b, c, alpha, beta, gamma
    :type cell:  list
    :param probability: thermal probability of the ellipsoid
    :type probability: float or int
    :param longest: not always the length is important. make to False to
                    get all three coordiantes of the ellipsoid axes.
    :type longest: boolean

    """
    from misc import A
    probability += 1
    # Uij is symmetric:
    if len(uvals) != 6:
        raise Exception('6 Uij values have to be supplied!')
    if len(cell) != 6:
        raise Exception('cell needs six parameters!')
    # orthogonalization matrix that transforms the fractional coordinates
    # with respect to a crystallographic basis system to coordinates
    # with respect to a Cartesian basis:
    A = A(cell).orthogonal_matrix
    Ucart = ufrac_to_ucart(A, cell, uvals)
    # print(Ucart)
    # E => eigenvalues, Q => eigenvectors:
    E, Q = mpm.eig(Ucart)
    # calculate vectors of ellipsoid axes  
    try:
        sqrt(E[0])
        sqrt(E[1])
        sqrt(E[2])
    except ValueError:
        print('\nEllipsoid is non positive definite!\n')
        return False
    v1 = mpm.matrix([Q[0, 0], Q[1, 0], Q[2, 0]])
    v2 = mpm.matrix([Q[0, 1], Q[1, 1], Q[2, 1]])
    v3 = mpm.matrix([Q[0, 2], Q[1, 2], Q[2, 2]])
    v1i = v1 * (-1)
    v2i = v2 * (-1)
    v3i = v3 * (-1)
    # multiply probability (usually 50%)
    e1 = sqrt(E[0]) * probability
    e2 = sqrt(E[1]) * probability
    e3 = sqrt(E[2]) * probability
    # scale axis vectors to eigenvalues 
    v1, v2, v3, v1i, v2i, v3i = v1 * e1, v2 * e2, v3 * e3, v1i * e1, v2i * e2, v3i * e3
    # find out which vector is the longest:
    length = mpm.norm(v1)
    v = 0
    if mpm.norm(v2) > length:
        length = mpm.norm(v2)
        v = 1
    elif mpm.norm(v3) > length:
        length = mpm.norm(v3)
        v = 2
    # move vectors back to atomic position
    atom = A * mpm.matrix(coords)
    v1, v1i = v1 + atom, v1i + atom
    v2, v2i = v2 + atom, v2i + atom
    v3, v3i = v3 + atom, v3i + atom
    # go back into fractional coordinates:
    a1 = cart_to_frac(v1, cell)
    a2 = cart_to_frac(v2, cell)
    a3 = cart_to_frac(v3, cell)
    a1i = cart_to_frac(v1i, cell)
    a2i = cart_to_frac(v2i, cell)
    a3i = cart_to_frac(v3i, cell)
    allvec = [[a1, a2, a3], [a1i, a2i, a3i]]
    if longest:
        # only the longest vector
        return [allvec[0][v], allvec[1][v]]
    else:
        # all vectors:
        return allvec


def ufrac_to_ucart(A, cell, uvals):
    U11, U22, U33, U23, U13, U12 = uvals
    U21 = U12
    U32 = U23
    U31 = U13
    Uij = mpm.matrix([[U11, U12, U13], [U21, U22, U23], [U31, U32, U33]])
    a, b, c, alpha, beta, gamma = cell
    V = vol_unitcell(*cell)
    # calculate reciprocal lattice vectors:
    astar = (b * c * sin(radians(alpha))) / V
    bstar = (c * a * sin(radians(beta))) / V
    cstar = (a * b * sin(radians(gamma))) / V
    # matrix with the reciprocal lattice vectors:
    N = mpm.matrix([[astar, 0, 0],
                    [0, bstar, 0],
                    [0, 0, cstar]])
    # Finally transform Uij values from fractional to cartesian axis system:
    Ucart = A * N * Uij * N.T * A.T
    return Ucart


def almost_equal(a, b, places=3):
    """
    Returns True or False if the number a and b are are equal inside the
    decimal places "places".
    :param a: a real number
    :type a: int/float
    :param b: a real number
    :type b: int/float
    :param places: number of decimal places
    :type places: int

    >>> almost_equal(1.0001, 1.0005)
    True
    >>> almost_equal(1.1, 1.0005)
    False
    >>> almost_equal(2, 1)
    False
    """
    return round(abs(a - b), places) == 0


def file_len(fname):
    """
    Returns the number of lines in a text file.

    >>> file_len('p21c.hkl')
    42976
    """
    i = 0
    with open(fname) as f:
        for i, l in enumerate(f, 1):
            pass
    return i


def get_overlapped_chunks(ring, size):
    """
    returns a list of chunks of size 'size' which overlap with one field.
    If the last chunk is smaller than size, the last 'size' chunks are returned as last chunk.
    "size" has to be larger than 3 to get reasonable results.
    >>> l = [1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 'a', 'b', 'c', 'd', 'e']
    >>> get_overlapped_chunks(l, 4)
    [[1, 2, 3, 4], [4, 5, 6, 7], [0, 7, 8, 9], [0, 'a', 'b', 'c'], ['b', 'c', 'd', 'e']]
    >>> get_overlapped_chunks(l, 3)
    [['c', 'd', 'e'], ['c', 'd', 'e'], ['c', 'd', 'e'], ['c', 'd', 'e'], ['c', 'd', 'e']]
    >>> get_overlapped_chunks(l, 5)
    [[1, 2, 3, 4, 5], [4, 5, 6, 7, 8], [0, 7, 8, 9, 'a'], [0, 'a', 'b', 'c', 'd'], ['a', 'b', 'c', 'd', 'e']]
    """
    chunks = []
    for i in range(0, len(ring) - size + 3, 3):
        chunk = ring[i:i + size]
        if len(chunk) < 4:
            chunk = ring[-size:]
        chunks.append(sorted(chunk))
    return chunks


def chunks(l, n):
    """
    returns successive n-sized chunks from l.

    >>> l = [1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 'a', 'b', 'c', 'd', 'e', 'f']
    >>> chunks(l, 5)
    [[1, 2, 3, 4, 5], [6, 7, 8, 9, 0], ['a', 'b', 'c', 'd', 'e'], ['f']]
    >>> chunks(l, 1)
    [[1], [2], [3], [4], [5], [6], [7], [8], [9], [0], ['a'], ['b'], ['c'], ['d'], ['e'], ['f']]
    >>> chunks(l, 50)
    [[1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 'a', 'b', 'c', 'd', 'e', 'f']]
    """
    out = []
    for i in range(0, len(l), n):
        out.append(l[i:i + n])
    return out


if __name__ == '__main__':
    import sys
    import doctest

    failed, attempted = doctest.testmod()  # verbose=True)
    if failed == 0:
        print('passed all {} tests!'.format(attempted))
