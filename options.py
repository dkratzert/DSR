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
try:
    from argparse import RawTextHelpFormatter
except(ImportError):
    print('\nYour Python version is incompatible with DSR! \nPlease use Python 2.7 or above.\n')
    sys.exit()
from argparse import ArgumentParser, SUPPRESS


__metaclass__ = type  # use new-style classes



class OptionsParser():
    '''
    This class uses the ArgumentParser module to parse the command line options.
    '''
    def __init__(self, program_name=''):
        self.progname = program_name
        self._options = self.parse_commandline()
        return

    def error(self):
        #print("\nPlease give one of the options as argument!\n")
        self.parser.print_help()
        sys.exit()

    @property
    def res_file(self):
        return self._options.res_file

    @property
    def external_restr(self):
        return self._options.external_restr

    @property
    def no_refine(self):
        return self._options.no_refine

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
    def import_grade(self):
        return self._options.import_grade

    @property
    def invert(self):
        return self._options.invert

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
    def all_options(self):
        return self._options



    def parse_commandline(self):
        '''parses the command line options and returns
           the command line options as dict'''
        # Options parser for the command line:
        sep_line = '\n------------------------------------------------'\
                    '--------------------------------\n'
        self.parser = ArgumentParser(prog='dsr', formatter_class=RawTextHelpFormatter,
        description='Disordered Structure Refinement (DSR)\n'
        +'\n'+self.progname
        +'\nExample DSR .res file command line:\n'
        +'\nREM DSR PUT/REPLACE "Fragment" WITH C1 C2 C3 ON Q1 Q2 Q3 PART 1 OCC -21 ='
        +'\n  RESI DFIX\n'
        +sep_line
        +'   PUT:     Just put the fragment source atoms here.\n'
        +'   REPLACE: Replace atoms of PART 0 in 1.2 A distance around target atoms.'
        +sep_line
        )
        self.parser.add_argument("-r", dest="res_file", metavar='"res file"', \
                                help="res file with DSR command", default=False)
        self.parser.add_argument("-re", dest="external_restr", metavar='"res file"', \
                                help="res file with DSR command (write restraints to external file)", default=False)
        self.parser.add_argument("-e", dest="export_fragment", metavar='"fragment"', \
                                help="export fragment as .res/.png file", default=False)
        self.parser.add_argument("-c", dest="export_clip", metavar='"fragment"', \
                                help="export fragment to clipboard", default=False)
        self.parser.add_argument("-t", dest="invert", action='store_true', \
                                help="inverts the current fragment", default=False)
        self.parser.add_argument("-i", dest="import_grade", metavar='"tgz file"', \
                                help="import a fragment from GRADE (needs .tgz file)", default=False)
        self.parser.add_argument("-ea", dest="export_all", action='store_true', \
                                help=SUPPRESS, default=False)
        self.parser.add_argument("-l", dest="list_db", action="store_true", \
                                help="list names of all database entries", default=False)
        self.parser.add_argument("-s", dest="search_string", metavar='"string"', \
                                help="search the database for a name", default=False)
        self.parser.add_argument("-n", dest="no_refine", action="store_true", \
                                help="do not refine after fragment transfer", default=False)
        return self.parser.parse_args()





if __name__ == '__main__':
    from dsr import program_name
    optparse = OptionsParser(program_name)
    print(optparse.res_file)
    print(optparse.parse_commandline())
    optparse.parser.print_help()


