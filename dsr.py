#/usr/bin/env python
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
import os
import logging
from dsrparse import DSR_Parser
from options import OptionsParser
from dbfile import global_DB, ImportGRADE, search_fragment_name,\
    print_search_results
from atomhandling import SfacTable, get_atomtypes, check_source_target,\
    set_final_db_sfac_types, replace_after_fit
from atomhandling import FindAtoms, NumberScheme
from resi import Resi
from restraints import ListFile, Lst_Deviations
from restraints import Restraints
from misc import remove_file
import misc
from afix import InsertAfix
from terminalsize import get_terminal_size
from refine import ShelxlRefine
import resfile
from cf3fit import CF3
from constants import width, sep_line

VERSION = '1.7.3'
# dont forget to change version in Innoscript file, spec file and deb file.

program_name = '\n'+((width/2)-9)*'-'+\
           ' D S R - v{} '.format(VERSION)+\
           ((width/2)-8)*'-'

# TODO and ideas:
'''
- try pure-python fft to make fourier maps and q-peaks.
- fit fragment without user interaction
'''



class DSR():
    '''
    main class
    '''
    def __init__(self, res_file_name=None, external_restr=None,
                 export_fragment=None, search_string=None,
                 export_clip=None, import_grade=None, export_all=None,
                 list_db=None, no_refine=None, invert=None):
        '''
        :param res_file_name:  name of the SHELXL res file, like 'p21c.res'
        :type res_file_name:   string
        :param external_restr: turn external restraints file on or off
        :type external_restr:  boolean
        :param export_fragment: name of fragment to be exported, like 'pfanion'
        :type export_fragment: string
        :param search_string:  search for this string in the database
        :type search_string:   string
        :param export_clip:    export fragment to clipboard, like 'pfanion'
        :type export_clip:     string
        :param import_grade:   import grade file to user database
        :type import_grade:    string
        :param export_all:     export all fragments at once
        :type export_all:      boolean
        :param list_db:        list database entries
        :type list_db:         boolean
        :param no_refine:      turn refinement off
        :type no_refine:       boolean
        :param invert:         invert the fragment during import, export or LS-fit
        :type invert:          boolean
        '''
        print(program_name)
        # options from the commandline options parser:
        self.options = OptionsParser()
        # vars() retrieves the options as dict, values() the values and any()
        # decides if any option is set.
        self.external = False
        if not res_file_name:
            self.res_file = self.options.res_file
        else:
            self.res_file = res_file_name
        if self.options.external_restr:
            self.external = True
            self.res_file = self.options.external_restr
        if external_restr:
            self.external = True
            self.res_file = external_restr
        if not export_fragment:
            self.export_fragment = self.options.export_fragment
        else:
            self.export_fragment = export_fragment
        if not export_clip:
            self.export_clip = self.options.export_clip
        else:
            self.export_clip = export_clip
        if not import_grade:
            self.import_grade = self.options.import_grade
        else:
            self.import_grade = import_grade
        if not export_all:
            self.export_all = self.options.export_all
        else:
            self.export_all = self.export_all
        if not list_db:
            self.list_db = self.options.list_db
        else:
            self.list_db = list_db
        if not no_refine:
            self.no_refine = self.options.no_refine
        else:
            self.no_refine = no_refine
        if not invert:
            self.invert = self.options.invert
        else:
            self.invert = invert
        if not search_string:
            self.search_string = self.options.search_string
        else:
            self.search_string = search_string
        if self.invert:
            if not any([self.res_file, self.external, self.import_grade,
                       self.export_clip, self.export_all, self.export_fragment]):
                self.options.error()
        #  List of database Fragments:
        if self.list_db:
            self.list_dbentries()
        if self.search_string:
            result = search_fragment_name(self.search_string)
            print_search_results(result)
            sys.exit()
        ## Export !all! fragments
        if self.export_all:
            self.export_all_fragments()
        ## Export one fragment
        if self.export_fragment:
            try:
                self.do_export_fragment()
            except() as e:
                print(e)
        if self.export_clip:
            try:
                self.export_to_clipboard()
            except() as e:
                print(e)
        ## Import a GRADE fragment
        if self.import_grade:
            self.import_from_grade()
        if not any(vars(self.options.all_options).values()+[self.res_file]):
            self.options.error()
        import time
        time1 = time.clock()
        self.main()
        time2 = time.clock()
        runtime = (time2 - time1)
        print('Runtime: {:>.1f} s'.format(runtime))
        print('\nDSR run complete.')

###############################################################################

    def do_export_fragment(self):
        '''
        Exports the current fragment.
        '''
        from export import Export
        gdb = global_DB(self.invert)
        export = Export(self.export_fragment, gdb, self.invert)
        export.write_res_file()
        sys.exit()

    def export_all_fragments(self):
        '''
        export all database entries at once
        '''
        from export import Export
        gdb = global_DB(self.invert)
        db = gdb.build_db_dict()
        dbnames = list(db.keys())
        for name in dbnames:
            export = Export(name, gdb, self.invert, self.export_all)
            export.write_res_file()
        sys.exit(1)

    def export_to_clipboard(self):
        '''
        Exports the current fragment to the clipboard.
        '''
        from export import Export
        gdb = global_DB(self.invert)
        export = Export(self.export_clip, gdb)
        export.export_to_clip()
        sys.exit(True)

    def list_dbentries(self):
        '''
        list all entries in the db.
        '''
        try:
            (width, height) = get_terminal_size()  # @UnusedVariable
        except():
            width = 80
        gdb = global_DB()
        db = gdb.build_db_dict()
        print('\n Entries found in the databases:\n')
        print(' Fragment         | Line | DB Name    | Full name, Comments ')
        print(sep_line)
        frags = sorted(db.keys())
        names_list = []
        num = 0
        for num, i in enumerate(frags):
            fragname = gdb.get_comment_from_fragment(i)
            names_list.append([i, fragname])
            line = ' {:<17}| {:<5}| {:<11}| {}'.format(
                    i, gdb.get_line_number_from_fragment(i),
                    gdb.get_db_name_from_fragment(i), fragname)
            print(line[:width - 1])
        try:
            if os.environ["DSR_DB_DIR"]:
                dbdir = os.environ["DSR_DB_DIR"]
        except(KeyError):
            dbdir = '.'
        print('\n {} Fragments in the database(s).'.format(num),
              '\n Feel free to add more fragments to "{}dsr_user_db.txt"' \
              '\n or mail them to dkratzert@gmx.de.'.format(dbdir + os.path.sep))
        for fragment in list(db.keys()):
            gdb.check_consistency(db[fragment], fragment)
            gdb.check_db_atom_consistency(db[fragment]['atoms'], fragment)
            gdb.check_db_header_consistency(db[fragment]['head'], db[fragment]['atoms'], fragment)
        sys.exit()

    def import_from_grade(self):
        '''
        imports a fragment from the GRADE webserver.
        '''
        mog = ImportGRADE(self.import_grade, self.invert)
        mog.write_user_database()
        sys.exit(1)


    def set_post_refine_cycles(self, shx, cycles):
        '''
        set the number of refinement cycles
        cycles must be a string
        '''
        try:
            shx.set_refinement_cycles(cycles)
        except(IndexError):
            print('Unable to set refinement cycles')

    def go_refine(self, shx):
        '''
        actually starts the fragment fit
        '''
        try:
            shx.run_shelxl()
        except() as e:
            print(e)
            sys.exit()

    def main(self):
        '''
        main object to run DSR as command line program
        '''
        #print(program_name)
        # The database content:
        basefilename = resfile.filename_wo_ending(self.res_file)
        gdb = global_DB(self.invert)
        rl = resfile.ResList(self.res_file)
        reslist = rl.get_res_list()
        find_atoms = FindAtoms(reslist)
        rle = resfile.ResListEdit(reslist, find_atoms)
        dsrp = DSR_Parser(reslist, rle)
        dsr_dict = dsrp.get_dsr_dict
        fvarlines = rle.find_fvarlines()
        fragment = dsrp.fragment.lower()
        if fragment in ['cf3', 'cf6', 'cf9']:
            dbhead = 'RESI CF3'
            db_residue_string = 'CF3'
            db_atom_types = ['C', 'F']
        else:
            dbhead = gdb.get_head_from_fragment(fragment)        # this is only executed once
            db_residue_string = gdb.get_resi_from_fragment(fragment)
            dbatoms = gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
            # the atomtypes of the dbentry as list e.g. ['C', 'N', ...]
            db_atom_types = get_atomtypes(dbatoms)
        sf = SfacTable(reslist, db_atom_types)
        sfac_table = sf.set_sfac_table()                 # from now on this sfac table is set
        resi = Resi(reslist, dsr_dict, dbhead, db_residue_string, find_atoms)
        # line where the dsr command is found in the resfile:
        dsr_line_number = dsrp.find_dsr_command(line=False)
        if fragment in ['cf3', 'cf6', 'cf9']:
            cf3 = CF3(rle, find_atoms, reslist, fragment, sfac_table,
                      basefilename, dsr_dict, resi, self.res_file)
            if fragment == 'cf3':
                cf3.cf3()
            if fragment == 'cf6':
                cf3.cf3('120')
            if fragment == 'cf9':
                cf3.cf9()
            print('\nFinised...')
            sys.exit()
        if dsrp.occupancy:
            rle.set_free_variables(dsrp.occupancy)
        fragline = gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
        dbhead = resi.remove_resi(dbhead)
        ### corrects the atom type according to the previous defined global sfac table:
        dbatoms = set_final_db_sfac_types(db_atom_types, dbatoms, sfac_table)

        ## Insert FRAG ... FEND entry:
        rle.insert_frag_fend_entry(dbatoms, fragline, fvarlines)

        print('Inserting {} into res File.'.format(fragment))
        if self.invert:
            print('Fragment inverted.')
        print('Source atoms:', ', '.join(dsrp.source))
        print('Target atoms:', ', '.join(dsrp.target))

        # several checks if the atoms in the dsr command line are consistent
        check_source_target(dsrp.source, dsrp.target, dbatoms)
        num = NumberScheme(reslist, dbatoms, resi.get_resinumber)
        # returns also the atom names if residue is active
        fragment_numberscheme = num.get_fragment_number_scheme()
        print('Fragment atom names:', ', '.join(fragment_numberscheme))
        dfix_head = ''
        if dsrp.dfix_active:
            restr = Restraints(fragment, gdb)
            dfix_12 = restr.get_formated_12_dfixes()
            dfix_13 = restr.get_formated_13_dfixes()
            flats = restr.get_formated_flats()
            dfix_head = dfix_12+dfix_13+flats
        afix = InsertAfix(reslist, dbatoms, db_atom_types, dbhead, dsr_dict,
                          sfac_table, find_atoms, fragment_numberscheme, dfix_head)
        afix_entry = afix.build_afix_entry(self.external, basefilename+'.dfix',
                                           resi)
        if dsr_line_number < fvarlines[-1]:
            print('\nWarning! The DSR command line MUST not appear before FVAR or the first atom in the .res file!')
            print('Can not proceed...\n')
            sys.exit()
        reslist[dsr_line_number] = reslist[dsr_line_number]+'\n'+afix_entry
        #reslist.insert(dsr_line_number+1, afix_entry)
        # write to file:
        shx = ShelxlRefine(reslist, basefilename, find_atoms)
        acta_lines = shx.remove_acta_card()
        shx.set_refinement_cycles('0')
        rl.write_resfile(reslist, '.ins')
        if self.no_refine:
            print('\nPerforming no fragment fit. Just prepared the .ins file for you.')
            return
        #  Refine with L.S. 0 to insert the fragment
        self.go_refine(shx)
        # Display the results from the list file:
        lf = ListFile(basefilename)
        lst_file = lf.read_lst_file()
        shx.check_refinement_results(lst_file)
        lfd = Lst_Deviations(lst_file)
        lfd.print_LS_fit_deviations()
        cell = rle.get_cell()
        # open res file again to restore 8 refinement cycles:
        rl = resfile.ResList(self.res_file)
        reslist = rl.get_res_list()
        if dsrp.command == 'REPLACE':
            reslist, find_atoms = replace_after_fit(rl, reslist, resi,
                                    fragment_numberscheme, cell)
        shx = ShelxlRefine(reslist, basefilename, find_atoms)
        shx.restore_acta_card(acta_lines)
        self.set_post_refine_cycles(shx, '8')
        shx.remove_afix()   # removes the afix 9
        # final resfile write:
        rl.write_resfile(reslist, '.res')


if __name__ == '__main__':
    '''main function'''
    #dsr = DSR(list_db=True)
    #dsr = DSR(export_fragment='toluene')
    #dsr = DSR(res_file_name='p21n_cf3.res')
    #sys.exit()
    #import cProfile
    try:
        #cProfile.run('dsr = DSR(res_file_name="p21c.res")', 'foo.profile')
        remove_file(misc.reportlog)
        dsr = DSR()
    except Exception as e:
        import platform
        logging.basicConfig(filename=misc.reportlog, filemode='w', level=logging.DEBUG)
        remove_file(misc.reportlog)
        logging.info('DSR version: {}'.format(VERSION))
        logging.info('Python version: {}'.format(sys.version))
        try:
            logging.info('Platform: {} {}, {}'.format(platform.system(),
                               platform.release(), ' '.join(platform.uname())))
        except:
            pass
        logger = logging.getLogger('dsr')
        ch = logging.StreamHandler()
        logger.addHandler(ch)
        print('\n\n')
        print('Congratulations! You found a bug in DSR. Please send the file\n'\
              ' "report-bug.log" and the .res file (if possible) to dkratzert@gmx.de\n')
        logger.exception(e)



