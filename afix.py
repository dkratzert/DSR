# -*- encoding: utf-8 -*-
# mÃ¶p
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <dkratzert@gmx.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
from __future__ import print_function

import fnmatch
import os
import sys

import misc
from constants import RESTRAINT_CARDS

__metaclass__ = type  # use new-style classes


def add_residue_to_dfix(dfix_head, resinum):
    """
    Add a residue to a list of DFIX/DANG restraints
    DFIX 1.234 C1 C2 -> DFIX 1.234 C1_4 C2_4
    >>> add_residue_to_dfix(['DFIX 1.456 C1 C2', 'DFIX 1.212 C3 C4'], 4)
    ['DFIX  1.456  C1_4  C2_4\\n', 'DFIX  1.212  C3_4  C4_4\\n']
    >>> add_residue_to_dfix(['DFIX 1.456 C1 C2', 'DFIX 1.212 C3 C4'], '5')
    ['DFIX  1.456  C1_5  C2_5\\n', 'DFIX  1.212  C3_5  C4_5\\n']
    >>> add_residue_to_dfix(['FLAT C6 C1 C2 C3', 'FLAT  C7  C1  C2  C3'], '2')
    ['FLAT  C6_2  C1_2  C2_2  C3_2\\n', 'FLAT  C7_2  C1_2  C2_2  C3_2\\n']
    """
    newhead = []
    for line in dfix_head:
        line = line.split()
        first = line[0]
        for num, item in enumerate(line):
            try:
                int(item[0])
            except:
                line[num] = line[num] + '_' + str(resinum)
                continue
        line = first + '  ' + '  '.join(line[1:]) + '\n'
        newhead.append(line)
    return newhead


def collect_all_restraints(reslist):
    """
    :return all_restraints: list
    collects all restraints in the resfile and returns a list with them
    [['RIGU_CF3', 'O1', '>', 'F9'], '...']
    """
    all_restraints = []
    for n, resline in enumerate(reslist):
        resline = resline.strip(' \n\r')
        # resline = resline.split()
        try:
            resline[:4]
        except:
            continue
        if resline[:4] in RESTRAINT_CARDS:
            # see for the next  lines if the lines continues with "=":
            line = 0
            while resline[-1] == '=':
                resline = resline[:-1] + reslist[n + line + 1]
                line += 1
                if not resline[-1] == '=':
                    break
                if line > 500:
                    break
            all_restraints.append(resline)
    return all_restraints


def remove_duplicate_restraints(reslist, restraints, residue_class=''):
    """
    removes restraints from the header which are already
    in the res-file.

    :param restraints:         database header (list of strings)
    :param residue_class:  SHELXL residue class
    :type residue_class:   string
    :param all_restraints: all restraints in the res file
    :type all_restraints:  list
    :return new_restr: list
    >>> dbhead = ["SADI 0.02 C1 C2 C2 C3 C3 C4", "SADI 0.04 C1 C3 C3 C5", "DFIX 1.45 C1 C2"]
    >>> all_restraints = ["SADI C1 C2 C2 C3 C3 C4", "SADI C1 C3 C3 C5", "DFIX C1 C2", "SADI C4 C5 C5 C6"]
    >>> remove_duplicate_restraints(all_restraints, dbhead)
    <BLANKLINE>
    Already existing restraints were not applied again.
    []

    >>> all_restraints = ["SADI 0.02 C1 C2 C2 C3 C3 C4", "SADI 0.04 C1 C3 C3 C5", "DFIX 1.45 C1 C2"]
    >>> dbhead = ["SADI C1 C2 C2 C3 C3 C4", "SADI C1 C3 C3 C5", "DFIX C1 C2", "SADI C4 C5 C5 C6"]
    >>> remove_duplicate_restraints(all_restraints, dbhead)
    <BLANKLINE>
    Already existing restraints were not applied again.
    ['SADI C4 C5 C5 C6']
    """
    modified = False
    new_restr = restraints[:]
    for num, line in enumerate(restraints):
        line = remove_stddev_from_restraint(line.split())
        for restr in collect_all_restraints(reslist):
            restr = restr.split()
            restr = remove_stddev_from_restraint(restr)
            if line == restr:
                new_restr[num] = ''
                modified = True
                break
    if modified:
        if residue_class:
            print('\nAlready existing restraints for residue "{}" were not '
                  'applied again.'.format(residue_class))
        else:
            print('\nAlready existing restraints were not applied again.')
    new_restr = [x for x in new_restr if x]
    return new_restr


def remove_stddev_from_restraint(restr):
    # type: (list) -> list
    """
    Parameters
    ----------
    restr: list of restraints

    Returns list of restraints without standars deviation
    -------

    >>> r = ['SADI', '0.02', 'C1', 'C2', 'C3', 'C4']
    >>> r2 = ['SADI', 'C1', 'C2', 'C3', 'C4']
    >>> remove_stddev_from_restraint(r)
    ['SADI', 'C1', 'C2', 'C3', 'C4']
    >>> remove_stddev_from_restraint(r2)
    ['SADI', 'C1', 'C2', 'C3', 'C4']
    """
    new = []
    # find out where the atoms begin (leave out numbers):
    for num, i in enumerate(restr):
        if i[0].isalpha():
            new.append(i)
    return new


def write_dbhead_to_file(dsrp, filename, dbhead, resi_class, resi_number):
    """
    write the restraints to an external file
    :param filename:     filename of database file
    :param dbhead:       database header
    :param resi_class:   SHELXL residue class
    :param resi_number:  SHELXL residue number
    :return filename:    full file name where restraints will be written
    """
    number = '1'
    files = []
    # find a unique number for the restraint file:
    for filen in misc.sortedlistdir('.'):
        if fnmatch.fnmatch(filen, 'dsr_*_' + filename):
            filenum = filen.split('_')
            if str.isdigit(filenum[1]):
                files.append(filenum[1])
    filepath, filename = os.path.split(os.path.abspath(filename))
    try:
        number = str(int(files[-1]) + 1)
    except IndexError:
        pass
    if not dsrp.resiflag:  # no residues
        filename = 'dsr_' + number + '_' + filename
    if dsrp.resiflag and resi_number:  # only residue number known
        filename = 'dsr_' + resi_class + '_' + resi_number + '_' + filename
    if dsrp.resiflag and not resi_number and resi_class:  # only residue class known
        filename = 'dsr_' + resi_class + '_' + filename
    if os.path.isfile(os.path.abspath(filename)):
        print('Previous restraint file found. Using restraints from "{}"'.format(filename))
        return filename
    try:
        dfix_file = open(os.path.join(filepath, filename), 'w')  # open the ins file
    except IOError:
        print('*** Unable to write restraints file! Check directory write permissions. ***')
        sys.exit(False)
    print('Restraints were written to "{}"'.format(os.path.join(filepath, filename)))
    for i in dbhead:  # modified reslist
        dfix_file.write("%s" % i)  # write the new file
    dfix_file.close()
    return filename


if __name__ == '__main__':
    import doctest

    failed, attempted = doctest.testmod()  # verbose=True)
    if failed == 0:
        print('passed all {} tests!'.format(attempted))
