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
import sys
import re
import os

import misc
from constants import sep_line


try:
    from argparse import RawTextHelpFormatter
except ImportError:
    print('\nYour Python version is incompatible with DSR! \nPlease use Python 2.7 or above.\n')
    sys.exit()
from argparse import ArgumentParser, SUPPRESS


__metaclass__ = type  # use new-style classes



class OptionsParser():
    """
    This class uses the ArgumentParser module to parse the command line options.
    """
    def __init__(self, versionline):
        self.versionline = versionline
        self.parser = ArgumentParser(prog='dsr', formatter_class=RawTextHelpFormatter,
                                     description='{}\nDisordered Structure Refinement (DSR)\n'.format(self.versionline)
                                     + '\nExample DSR .res file command line:'
                                     + '\nREM DSR PUT/REPLACE "Fragment" WITH C1 C2 C3 ON Q1 Q2 Q3 PART 1 OCC -21 ='
                                     + '\n  RESI DFIX\n'
                                     + sep_line+'\n'
                                     + '   PUT:     Just put the fragment source atoms here.\n'
                                     + '   REPLACE: Replace atoms of PART 0 in 1.3 A distance around target atoms.\n'
                                     + sep_line
                                     )
        self.parser.add_argument("-r", dest="res_file", metavar='"res file"', nargs='+',
                                 help="res file with DSR command", default=False)
        self.parser.add_argument("-re", dest="external_restr", metavar='"res file"', nargs='+',
                                 help="res file with DSR command (write restraints to external file)", default=False)
        self.parser.add_argument("-e", dest="export_fragment", metavar='"fragment"',
                                 help="export fragment as .res/.png file", default=False)
        self.parser.add_argument("-c", dest="export_clip", metavar='"fragment"',
                                 help="export fragment to clipboard", default=False)
        self.parser.add_argument("-t", dest="invert", action='store_true',
                                 help="inverts the current fragment", default=False)
        self.parser.add_argument("-i", dest="import_grade", metavar='"tgz file"',
                                 help="import a fragment from GRADE (needs .tgz file)", default=False)
        self.parser.add_argument("-l", dest="list_db", action="store_true",
                                 help="list names of all database entries", default=False)
        self.parser.add_argument("-s", dest="search_string", metavar='"string"',
                                 help="search the database for a name", default=False)
        self.parser.add_argument("-g", dest="rigid_group", help="keep group rigid (no restraints)",
                                 action="store_true", default=False)
        self.parser.add_argument("-u", dest="selfupdate", help="Update DSR to the most current version",
                                 action="store_true", default=False)
        self.parser.add_argument("-ea", dest="export_all", action='store_true',
                                 help=SUPPRESS, default=False)
        self.parser.add_argument("-lc", dest="list_db_csv", action='store_true',
                                 help=SUPPRESS, default=False)
        self.parser.add_argument("-x", dest="search_extern",
                                 help=SUPPRESS, default=False)
        self.parser.add_argument("-ah", dest="head_for_gui",
                                 help=SUPPRESS, default=False)
        # with nargs='+' it accepts space in path and returns a list:
        self.parser.add_argument("-shx", dest="shelxl_ex", nargs='+',
                                 help=SUPPRESS, default=False)
        self.parser.add_argument("-n", dest="no_refine", action="store_true",
                                 help="do not refine after fragment transfer", default=False)
        self._options = self.parser.parse_args()

    def error(self):
        self.parser.print_help()
        sys.exit()

    @property
    def res_file(self):
        if self._options.res_file:
            rpath = r' '.join(self._options.res_file)
        else:
            return False
        try:
            rpath = os.path.normpath(rpath)
        except:
            rpath = None
        return rpath

    @property
    def external_restr(self):
        if self._options.external_restr:
            erpath = r' '.join(self._options.external_restr)
        else:
            return False
        try:
            erpath = os.path.normpath(erpath)
        except:
            erpath = None
        return erpath

    @property
    def no_refine(self):
        return self._options.no_refine

    @property
    def rigid_group(self):
        return self._options.rigid_group

    @property
    def export_fragment(self):
        return self._options.export_fragment

    @property
    def export_clip(self):
        return self._options.export_clip

    @property
    def export_all(self):
        return self._options.export_all

    @property
    def list_db(self):
        return self._options.list_db

    @property
    def list_db_csv(self):
        return self._options.list_db_csv

    @property
    def import_grade(self):
        return self._options.import_grade

    @property
    def selfupdate(self):
        return self._options.selfupdate

    @property
    def head_for_gui(self):
        frag = False
        if self._options.head_for_gui:
            frag = self._options.head_for_gui.lower()
        return frag

    @property
    def invert(self):
        return self._options.invert

    @property
    def shelxl_ex(self):
        '''
        Option to define the path of the shelxl executable
        '''
        spath = ''
        if self._options.shelxl_ex:
            spath = r' '.join(self._options.shelxl_ex)
        else:
            return False
        if not os.access(spath, os.X_OK):
            try:
                spath = misc.which(spath)
                return spath[0]
            except IndexError:
                return ''
        return spath

    @property
    def search_string(self):
        if not self._options.search_string:
            return None
        # search characters allowed: a-z A-Z 0-9 _ - , () {} [] ' " + * | = .
        chars = re.match(r'^[\w\-,\(\)\[\]\{\}\'\"\+\*\|\=\.]+$', 
                         self._options.search_string)
        if not chars:
            print('Characters not allowed for searching.')
            sys.exit()
        else:
            return self._options.search_string

    @property
    def search_extern(self):
        if not self._options.search_extern:
            return None
        # search characters allowed: a-z A-Z 0-9 _ - , () {} [] ' " + * | = .
        chars = re.match(r'^[\w\-,\(\)\[\]\{\}\'\"\+\*\|\=\.]+$', 
                         self._options.search_extern)
        if not chars:
            print('Characters not allowed for searching.')
            sys.exit()
        else:
            return self._options.search_extern

    @property
    def all_options(self):
        return self._options


if __name__ == '__main__':
    optparse = OptionsParser('188')
    print(optparse.res_file)
    #print(optparse.parse_commandline())
    optparse.parser.print_help()


