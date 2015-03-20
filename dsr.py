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
from dbfile import global_DB, ImportGRADE
from resfile import ResList, ResListEdit, filename_wo_ending
from atomhandling import SfacTable, get_atomtypes, check_source_target,\
    set_final_db_sfac_types
from atomhandling import FindAtoms, NumberScheme
from resi import Resi
from restraints import ListFile, Lst_Deviations
from restraints import Restraints
from misc import remove_file
import misc
from afix import InsertAfix
from terminalsize import get_terminal_size
from refine import ShelxlRefine

VERSION = '1.6.0'
# dont forget to change version in Innoscript file, spec file and deb file.
program_name = '\n-----------------------------'\
           ' D S R - v{}' \
           ' ----------------------------------'.format(VERSION)

# TODO and ideas:
# - In replace mode, check for atoms in PART 0 which are near the fitting fragment
#   e.g. below 1.4 Angstroms and delete them.
#   To prevent deletion of needed atoms, only delete atoms that are shorter
#   away than the covalence distance of the atom pair. maybe 2/3 of the distance?
# -detect collinear atoms
# -import also from pdb, dfix, obprop alone.
# -debian package: /usr/src/packages/BUILD # dpkg-deb --build dsr



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
        # options from the commandline options parser:
        self.options = OptionsParser(program_name)
        # vars() retrieves the options as dict, values() the values and any()
        # decides if any option is set.
        if not any(vars(self.options.all_options).values()):
            self.options.error()
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
                       self.export_clip, self.export_all]):
                self.options.error()
        #  List of database Fragments:
        if self.list_db:
            self.list_dbentries()
        if self.search_string:
            result = self.search_fragment_name()
            self.print_search_results(result)
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
                self.export_to_clip()
            except() as e:
                print(e)
        ## Import a GRADE fragment
        if self.import_grade:
            self.import_from_grade()

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

    def export_to_clip(self):
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
        print(' ----------------------------------------'
                '-----------------------------------')
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


    def search_fragment_name(self):
        '''
        searches the Name: comments in the database for a given name
        '''
        from misc import dice_coefficient
        #from misc import levenshtein
        gdb = global_DB()
        db = gdb.build_db_dict()
        frags = sorted(db.keys())
        names_list = []
        for i in frags:
            fragname = gdb.get_comment_from_fragment(i)
            line_number = gdb.get_line_number_from_fragment(i)
            names_list.append([i, fragname, line_number])
        search_results = {}
        for i in names_list:
            db_entry = i[1]
            #Levenshtein gibt bei kurzen Suchstrings zu schlechte Ergebnisse:
            #coefficient = levenshtein(self.search_string, db_entry)
            coefficient = dice_coefficient(self.search_string, db_entry)
            search_results[coefficient] = i
        # select the best 4 results:
        selected_results = [search_results[i] for i in sorted(search_results)[0:4]]
        return selected_results


    def print_search_results(self, results):
        '''
        prints the results of a database search to screen and exit.
        results are
        '''
        print('\n\n Found following database entries:\n')
        print(' Fragment          | Full name, Comments                      | Line number')
        print(' ---------------------------------------------------------------------------')
        for line in results:
            print(' {:15s}   | {:40s} | {}'.format(line[0], line[1], line[2]))
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


    def replacemode(self, res_target_atoms, rle, reslist, sfac_table):
        '''
        Target atoms are being replaced if this is executed
        '''
        for i in res_target_atoms:
            if '_' in i:
                print('\nDo you really want to REPLACE atom {} inside a residue?'.format(i))
                print('This will very likely damage something.\n')
                break
        fa = FindAtoms(reslist)
        print('Replace mode active.')
        target_lines = fa.get_atom_line_numbers(res_target_atoms)
        for i in target_lines:
            i = int(i)
            rle.remove_line(i, rem=False, remove=False, frontspace=True)
        h_delcount = fa.remove_adjacent_hydrogens(res_target_atoms, sfac_table)
        if h_delcount:
            return target_lines+h_delcount
        else:
            return target_lines


    def go_refine(self, shx):
        '''
        actually starts the fragment fit
        '''
        try:
            shx.run_shelxl()
        except() as e:
            print(e)
            sys.exit()


    def prepare_no_refine(self, shx, rl, reslist):
        shx.remove_acta_card()
        shx.set_refinement_cycles('0')
        rl.write_resfile(reslist, '.ins')
        print('\nPerforming no fragment fit. Just prepared the .ins file for you.')



    def replace_after_fit(self, rl, reslist, resi, fragment_numberscheme, cell):
        '''
        deletes the atoms in replace mode that are < 1.2 A near the fragment atoms
        '''
        find_atoms = FindAtoms(reslist)
        if resi.get_resinumber:
            frag_at = []
            for i in fragment_numberscheme:
                at = i + '_{}'.format(resi.get_resinumber)
                frag_at.append(at)
        else:
            frag_at = fragment_numberscheme
        atoms_to_delete = find_atoms.remove_near_atoms(frag_at, cell)
        print('Replacing following atoms (< 1.2 A near fragment):\n', 
              ' '.join(atoms_to_delete))
        target_lines = find_atoms.get_atom_line_numbers(atoms_to_delete)
        rle = ResListEdit(reslist, find_atoms)
        for i in target_lines:
            i = int(i)
            rle.remove_line(i, rem=False, remove=False, frontspace=True)
        
        rl.write_resfile(reslist, '.res')
        reslist = rl.get_res_list()
        return reslist, find_atoms

    def main(self):
        '''
        main object to run DSR as command line program
        '''
        print(program_name)
        # The database content:
        gdb = global_DB(self.invert)
        rl = ResList(self.res_file)
        reslist = rl.get_res_list()
        find_atoms = FindAtoms(reslist)
        rle = ResListEdit(reslist, find_atoms)
        dsrp = DSR_Parser(reslist, rle)
        dsr_dict = dsrp.parse_dsr_line()
        fvarlines = rle.find_fvarlines()
        if dsrp.occupancy:
            rle.set_free_variables(dsrp.occupancy, fvarlines)

        fragment = dsrp.fragment.lower()
        fragline = gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
        dbatoms = gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
        dbhead = gdb.get_head_from_fragment(fragment)        # this is only executed once
        db_residue_string = gdb.get_resi_from_fragment(fragment)
        db_atom_types = get_atomtypes(dbatoms)                 # the atomtypes of the dbentry as list e.g. ['C', 'N', ...]
        resi = Resi(reslist, dsr_dict, dbhead, db_residue_string, find_atoms)
        dbhead = resi.remove_resi(dbhead)
        sf = SfacTable(reslist, db_atom_types)
        sfac_table = sf.set_sfac_table()                 # from now on this sfac table is set

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
        basefilename = filename_wo_ending(self.res_file)
        num = NumberScheme(reslist, dbatoms, resi.get_resinumber)
        # returns also the atom names if residue is active
        fragment_numberscheme = num.get_fragment_number_scheme()
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
        # line where the dsr command is found in the resfile:
        dsr_line_number = dsrp.find_dsr_command(line=False)
        if dsr_line_number < fvarlines[-1]:
            print('\nWarning! The DSR command line MUST not appear before FVAR or the first atom in the .res file!')
            print('Can not proceed...\n')
            sys.exit()
        reslist[dsr_line_number] = reslist[dsr_line_number]+'\n'
        reslist.insert(dsr_line_number+1, afix_entry)

        ##### comment out all target atom lines in replace mode:
        if dsrp.command == 'REPLACE':
            self.replacemode(dsrp.target, rle, reslist, sfac_table)

        # write to file:
        shx = ShelxlRefine(reslist, basefilename, find_atoms)
        if self.no_refine:
            self.prepare_no_refine(shx, rl, reslist)
            sys.exit()
        shx.set_refinement_cycles('0')
        shx.remove_acta_card()
        rl.write_resfile(reslist, '.ins')
        #  Refine with L.S. 0 to insert the fragment
        self.go_refine(shx)
        # Display the results from the list file:
        lf = ListFile(basefilename)
        lst_file = lf.read_lst_file()
        lfd = Lst_Deviations(lst_file)
        lfd.print_LS_fit_deviations()
        cell = rle.get_cell()

        # open res file again to restore 8 refinement cycles:
        rl = ResList(self.res_file)
        reslist = rl.get_res_list()
        if dsrp.command == 'REPLACE':
            reslist, find_atoms = self.replace_after_fit(rl, reslist, resi, 
                                                    fragment_numberscheme, cell)
            
        shx = ShelxlRefine(reslist, basefilename, find_atoms)
#        shx.restore_acta_card()
        shx.check_refinement_results(lst_file)
        self.set_post_refine_cycles(shx, '8')
        shx.remove_afix()   # removes the afix 9

        # final resfile write:
        rl.write_resfile(reslist, '.res')


if __name__ == '__main__':
    '''main function'''
    #dsr = DSR(list_db=True)
    #dsr = DSR(res_file_name='p21c.res', external_restr=True)
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



