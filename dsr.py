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
from options import OptionsParser
from constants import width, sep_line
from misc import reportlog, remove_file, find_line,\
    remove_line
from dbfile import global_DB, search_fragment_name
from dsrparse import DSR_Parser
from dbfile import ImportGRADE, print_search_results
from atomhandling import SfacTable, get_atomtypes, check_source_target,\
    set_final_db_sfac_types, replace_after_fit
from atomhandling import FindAtoms, NumberScheme
from resi import Resi
from restraints import ListFile, Lst_Deviations
from restraints import Restraints
from afix import InsertAfix
from terminalsize import get_terminal_size
from refine import ShelxlRefine
import resfile
from cf3fit import CF3
from os.path import expanduser


VERSION = '187'
# dont forget to change version in Innoscript file, spec file and deb file.

program_name = '\n'+((width//2)-9)*'-'+\
                ' D S R - v{} '.format(VERSION)+\
                    ((width//2)-8)*'-'

# TODO and ideas:
'''

'''

class DSR():
    '''
    main class
    '''
    def __init__(self, res_file_name=None, external_restr=None, atom_coordinates=None,
                 export_fragment=None, search_string=None, search_extern=None,
                 export_clip=None, import_grade=None, export_all=None, rigid=None,
                 list_db=None, no_refine=None, invert=None, list_db_csv=None, head_csv=None):

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
        :param export_all:     hidden option, export all fragments at once
        :type export_all:      boolean
        :param list_db:        list database entries
        :type list_db:         boolean
        :param list_db_csv:    hidden option, list database entries in machine readyble form
        :type list_db_csv:     boolean
        :param no_refine:      turn refinement off
        :type no_refine:       boolean
        :param invert:         invert the fragment during import, export or LS-fit
        :type invert:          boolean
        '''
        # options from the commandline options parser:
        self.options = OptionsParser(program_name)
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
        if not list_db_csv:
            self.list_db_csv = self.options.list_db_csv
        else:
            self.list_db_csv = list_db_csv
        if not no_refine:
            self.no_refine = self.options.no_refine
        else:
            self.no_refine = no_refine
        if not invert:
            self.invert = self.options.invert
        else:
            self.invert = invert
        if not rigid:
            self.rigid = self.options.rigid_group
        else:
            self.rigid = rigid
        if not search_string:
            self.search_string = self.options.search_string
        else:
            self.search_string = search_string
        if not atom_coordinates:
            self.frag_for_gui = self.options.frag_for_gui
        if not search_extern:
            self.search_extern = self.options.search_extern
        else:
            self.search_extern = search_extern
        if not head_csv:
            self.head_csv = self.options.head_for_gui
        else:
            self.head_csv = head_csv
        if self.invert:
            if not any([self.res_file, self.external, self.import_grade,
                       self.export_clip, self.export_all, self.export_fragment]):
                self.options.error()
        if self.frag_for_gui:
            self.export_to_gui()
        if self.head_csv:
            self.head_to_gui()
        #  List of database Fragments:
        if self.list_db_csv:
            gdb = global_DB()
            frags = gdb.list_fragments()
            print('DSR version:', VERSION)
            for i in frags:
                print('{};;{};;{};;{}'.format(i[0], i[3], i[1], i[2]))
            sys.exit()  
        if self.search_extern:
            result = search_fragment_name(self.search_extern)
            for i in result:
                print('{};;{};;{};;{}'.format(i[0], i[1], i[2], i[3]))
            sys.exit()
        print(program_name)
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
        if not any(list(vars(self.options.all_options).values())+[self.res_file]):
            self.options.error()
        import time
        time1 = time.clock()
        self.main()
        time2 = time.clock()
        runtime = (time2 - time1)
        print('Runtime: {:>.1f} s'.format(runtime))
        print('DSR run complete.')

###############################################################################

    def head_to_gui(self):
        '''
        Exports current fragment header and atoms to the GUI
        '''
        from export import Export
        atoms = []
        helpmsg = "Please ask daniel.kratzert@ac.uni-freiburg.de for help."
        try:
            gdb = global_DB(self.invert, fragment=self.head_csv)
        except Exception as e:  # @UnusedVariable
            print("Initializing the database failed.")
            print(helpmsg)
            #print(e)
            sys.exit()
        try:
            export = Export(self.head_csv, gdb, self.invert)
        except:
            print("Unable to export informations from DSR.")
            sys.exit()
        try:
            atoms = export.export_to_gui()
        except:
            print("Could not get atom information.")
            print(helpmsg)
        print('<atoms>\n', atoms, '\n</atoms>')
        # prints most of the needed info:
        gdb.get_head_for_gui(self.head_csv)
        sys.exit()

    def export_to_gui(self):
        '''
        Exports the current fragment atoms to the GUI.
        '''
        from export import Export
        gdb = global_DB(self.invert)
        self.export_fragment = self.frag_for_gui
        export = Export(self.export_fragment, gdb, self.invert)
        atoms = export.export_to_gui()
        if not atoms:
            sys.exit()
        sys.exit()
    
    def do_export_fragment(self):
        '''
        Exports the current fragment to a res file.
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
        dbdir = expanduser('~')
        gdb = global_DB()
        fraglist = gdb.list_fragments()
        fragnames = []
        try:
            (width, height) = get_terminal_size()  # @UnusedVariable
        except():
            width = 80
        print('\n Entries found in the databases:\n')
        print(' Fragment         | Line | DB Name    | Full name, Comments ')
        print(sep_line)
        for num, line in enumerate(fraglist):
            fragnames.append(line[0])
            line = ' {:<17}| {:<5}| {:<11}| {}'.format(*line)
            print(line[:width - 1])
        print('\n {} Fragments in the database(s).'.format(num),
              '\n Feel free to add more fragments to "{}dsr_user_db.txt"' \
              '\n and please mail them to dkratzert@gmx.de.'.format(dbdir + os.path.sep))
        for fragment in fragnames:
            gdb.check_consistency(fragment)
            gdb.check_db_atom_consistency(fragment)
            gdb.check_db_header_consistency(fragment)
            gdb.check_sadi_consistence(fragment)
        sys.exit()

    def import_from_grade(self):
        '''
        imports a fragment from the GRADE webserver.
        '''
        mog = ImportGRADE(self.import_grade, self.invert)
        mog.write_user_database()
        sys.exit(1)

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
        if not basefilename:
            print('Illegal option')
            sys.exit()
        gdb = global_DB(self.invert)
        rl = resfile.ResList(self.res_file)
        reslist = rl.get_res_list()
        if len(reslist) == 0:
            print("The input file is empty. Can not proceed!")
            sys.exit()
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
                      basefilename, dsr_dict, resi, self.res_file, self.options)
            if fragment == 'cf3':
                cf3.cf3()
            if fragment == 'cf6':
                cf3.cf3('120')
            if fragment == 'cf9':
                cf3.cf9()
            print('\nFinished...')
            sys.exit()
        # checks have to be after CF3, CF6 etc.
        gdb.check_consistency(fragment)
        gdb.check_db_atom_consistency(fragment)
        gdb.check_db_header_consistency(fragment)
        gdb.check_sadi_consistence(fragment)
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
                          sfac_table, find_atoms, fragment_numberscheme, self.options, dfix_head)
        afix_entry = afix.build_afix_entry(self.external, basefilename+'.dfix', resi)
        if dsr_line_number < fvarlines[-1]:
            print('\nWarning! The DSR command line MUST NOT appear before FVAR or the first atom in the .res file!')
            print('Can not proceed...\n')
            sys.exit()
        reslist[dsr_line_number] = reslist[dsr_line_number]+'\n'+afix_entry
        #reslist.insert(dsr_line_number+1, afix_entry)
        # write to file:
        shx = ShelxlRefine(reslist, basefilename, find_atoms, self.options)
        acta_lines = shx.remove_acta_card()
        cycles = shx.get_refinement_cycles
        shx.set_refinement_cycles('0')  # also sets shx.cycles to current value
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
        # remove the "REM " instriction bevore the +dfixfile instruction
        plusline = find_line(reslist, "REM "+afix.rand_id_dfx)
        if plusline:
            reslist[plusline-1] = reslist[plusline-1][4:]
            remove_line(reslist, plusline, remove=True)
        if dsrp.command == 'REPLACE':
            reslist, find_atoms = replace_after_fit(rl, reslist, resi,
                                    fragment_numberscheme, cell)
        shx = ShelxlRefine(reslist, basefilename, find_atoms, self.options)
        shx.restore_acta_card(acta_lines)
        try:
            if cycles != None:
                shx.set_refinement_cycles(cycles) # restores last LS value
            else:
                shx.set_refinement_cycles(8)
        except(IndexError):
            print('Unable to set refinement cycles')
        if not self.options.rigid_group:
            shx.remove_afix(afix.rand_id_afix)   # removes the afix 9
        # final resfile write:
        rl.write_resfile(reslist, '.res')


if __name__ == '__main__':
    '''main function'''
    #dsr = DSR(list_db=True)
    #dsr = DSR(export_fragment='toluene')
    #dsr = DSR(res_file_name='p21n_cf3.res')
    #sys.exit()
    """
    import cProfile
    import pstats
    cp = cProfile.Profile()
    cp.enable(subcalls=True, builtins=True)
    """
    try:
        remove_file(reportlog)
        dsr = DSR()
    except Exception as e:
        import platform
        import logging
        logging.basicConfig(filename=reportlog, filemode='w', level=logging.DEBUG)
        remove_file(reportlog)
        logging.info('DSR version: {}'.format(VERSION))
        logging.info('Python version: {}'.format(sys.version))
        try:
            logging.info('Platform: {} {}, {}'.format(platform.system(),
                               platform.release(), ' '.join(platform.uname())))
        except:
            print("Can not write logfile")
            pass
        logger = logging.getLogger('dsr')
        ch = logging.StreamHandler()
        logger.addHandler(ch)
        print('\n')
        print('Congratulations! You found a bug in DSR. Please send the file\n'\
              ' "report-bug.log" and the .res file (if possible) to dkratzert@gmx.de\n'
              )
        #.format(os.path.dirname(os.path.realpath(reportlog))+os.sep ))
        logger.exception(e)
    """
    cp.disable()
    pstats.Stats(cp).sort_stats('cumtime').print_stats(30)
    pstats.Stats(cp).sort_stats('cumtime').print_callers(30)
    pstats.Stats(cp).sort_stats('cumtime').print_callees(30)
    """
    
