# /usr/bin/env python
# -*- encoding: utf-8 -*-
# möp
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <dkratzert@gmx.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
import re

from src.terminalsize import get_terminal_size

SHX_CARDS = ('TITL', 'CELL', 'ZERR', 'LATT', 'SYMM', 'SFAC', 'UNIT', 'LIST', 'L.S.', 'CGLS',
             'BOND', 'FMAP', 'PLAN', 'TEMP', 'ACTA', 'CONF', 'SIMU', 'RIGU', 'WGHT', 'FVAR',
             'DELU', 'SAME', 'DISP', 'LAUE', 'REM ', 'MORE', 'TIME', 'END', 'HKLF', 'OMIT',
             'SHEL', 'BASF', 'TWIN', 'EXTI', 'SWAT', 'HOPE', 'MERG', 'SPEC', 'RESI', 'MOVE',
             'ANIS', 'AFIX', 'HFIX', 'FRAG', 'FEND', 'EXYZ', 'EADP', 'EQIV', 'CONN', 'BIND',
             'FREE', 'DFIX', 'BUMP', 'SADI', 'CHIV', 'FLAT', 'DEFS', 'ISOR', 'NCSY', 'SUMP',
             'BLOC', 'DAMP', 'STIR', 'MPLA', 'RTAB', 'HTAB', 'SIZE', 'WPDB', 'GRID', 'MOLE',
             'XNPD', 'REST', 'CHAN', 'FLAP', 'RNUM', 'SOCC', 'PRIG', 'WIGL', 'RANG', 'TANG',
             'ADDA', 'STAG', 'NEUT', 'ABIN', 'ANSC', 'ANSR', 'NOTR', 'TWST', 'PART', 'DANG',
             'BEDE', 'LONE', 'REM', 'END')

RESTRAINT_CARDS = ('SIMU', 'RIGU', 'DELU', 'SAME', 'FREE', 'DFIX', 'BUMP', 'HFIX', 'BIND',
                   'SADI', 'CHIV', 'FLAT', 'DEFS', 'ISOR', 'NCSY', 'DANG', 'EADP', 'EXYZ')

# restraints regarding distances only:
DIST_RESTRAINT_CARDS = ('SAME', 'SADI', 'FREE', 'DFIX', 'BUMP', 'CHIV', 'FLAT', 'HFIX', 'DEFS',
                        'ISOR', 'NCSY', 'DANG')

# SHELXL atom definition:
# start at line beginning, one to 4 letters, zero to 3 digits,
# one or more whitespaces, one digit, one or more whitespaces, one or more digits,
# a dot, one or more digits, ...

atomregex = re.compile(
    r'^(([A-Za-z]{1,3})[0-9]{0,3}[A-Za-z\'\"]{0,2}\d{0,1})\s+[-]?[0-9]{1,2}\s+[-+]?[0-9]*\.?[0-9]+\s+[-+]?[0-9]*\.?[0-9]+\s+[-+]?[0-9]*\.?[0-9]+')

try:
    (width, height) = get_terminal_size()  # @UnusedVariable
except():
    width = 80
sep_line = (width - 1) * '-'

isoatomstr = '{:<5.5s} {:<3}{:>10.6f}  {:>10.6f}  {:>9.6f}  {:>9.5f}  {:>9.5f}'
