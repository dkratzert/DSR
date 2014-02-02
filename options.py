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
try:
    from argparse import RawTextHelpFormatter
except(ImportError):
    print('\nYour python version is incompatible with DSR! \nPlease use Python 2.7 or above.\n')
    sys.exit()
from argparse import ArgumentParser, SUPPRESS


__metaclass__ = type  # use new-style classes



class OptionsParser():
    '''
    This class uses the ArgumentParser module to parse the command line options.
    The options can be derived with:
    optparse = OptionsParser()
    optparse.res_file
    optparse.fragment
    '''
    def __init__(self, progname=''):
        self.progname = progname
        self._options = self.parse_commandline()
        # check if at least one option is given:
        if not self._options.res_file\
        and not self._options.export_fragment\
        and not self._options.list_db\
        and not self._options.export_all\
        and not self._options.import_grade\
        and not self._options.no_refine:
            self.error()
        
        #and not self._options.debug\
        
    def error(self):
        print("\nPlease give one of the options as argument!\n")
        self.parser.print_help()
        sys.exit()

    @property
    def res_file(self):
        return self._options.res_file
    
    @property
    def no_refine(self):
        return self._options.no_refine
    
    @property
    def export_fragment(self):
        return self._options.export_fragment
    
    @property
    def export_all(self):
        return self._options.export_all
    
    @property
    def debug(self):
        return self._options.debug
    
    @property
    def list_db(self):
        return self._options.list_db
    
    @property
    def import_grade(self):
        return self._options.import_grade
    
    
    def parse_commandline(self):
        '''parses the command line options and returns the command line options as dict'''
        # Options parser for the command line:
        sep_line = '\n--------------------------------------------------------------------------------\n'
        self.parser = ArgumentParser(prog='dsr', formatter_class=RawTextHelpFormatter, 
        description='Disordered solvent refinement (DSR)\n'
        +'\n'+self.progname
        +'\nExample DSR res-file command line:\n'
        +'\nREM DSR PUT/REPLACE "Fragment" WITH C1 C2 C3 ON Q1 Q2 Q3 PART 1 OCC 21.0 ='
        +'\n  RESI 1 BENZ\n'
        +sep_line
        +'   PUT:     Just put the fragment source atoms here.\n'
        +'   REPLACE: Replace existing target atoms or q-peaks.'
        +sep_line
        )
        self.parser.add_argument("-r", dest="res_file", metavar='"res file"', help="res file with DSR command", default=False)
        self.parser.add_argument("-e", dest="export_fragment", metavar='"fragment"', help="export fragment from the database", default=False)
        self.parser.add_argument("-i", dest="import_grade", metavar='"tgz file"', help="import a fragment from GRADE (needs .tgz file)", default=False)
        self.parser.add_argument("-ea", dest="export_all", action='store_true', help=SUPPRESS, default=False)
      #  self.parser.add_argument("-d", dest="debug", action="store_true", help="Enables printing of debug messages", default=False)
        self.parser.add_argument("-l", dest="list_db", action="store_true", help="list names of all database entries", default=False)
        self.parser.add_argument("-n", dest="no_refine", action="store_true", help="do not refine after fragment transfer", default=False)
        return self.parser.parse_args()
        




if __name__ == '__main__':
    from dsr import progname
    optparse = OptionsParser(progname)
    print(optparse.res_file)
    optparse.parser.print_help()

    
