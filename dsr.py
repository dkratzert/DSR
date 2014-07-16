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
from dsrparse import DSR_Parser
from options import OptionsParser
from dbfile import global_DB, ImportGRADE
from resfile import ResList, ResListEdit, filename_wo_ending
from atomhandling import SfacTable, get_atomtypes, check_source_target
from atomhandling import FindAtoms, NumberScheme, Elem_2_Sfac
from afix import InsertAfix, write_dbhead_to_file
from refine import ShelxlRefine
from resi import Resi
#from misc import find_line
from restraints import ListFile, Lst_Deviations, format_atom_names
from restraints import Restraints, Adjacency_Matrix

# TODO and ideas:
# -FLAT: see for every ring atom if there is also an attached atom in the same plane
# -To the manual: Avogradro/res-file -> rename -> mercury -> mol2-file -> GRADE
# -To manual: what to do if something does not work.
# -try to iterate all target atoms if fit is really bad
# -select atoms in shelxle > export to db
# -If second part of a disordered atom is very close to the first, make EADP?
# -detect empty residues and parts after atom deletion
# -import also from pdb, dfix, obprop alone.
# -debian package: /usr/src/packages/BUILD # dpkg-deb --build dsr


VERSION = '1.5.7'
progname = '\n-----------------------------'\
           ' D S R - v{}' \
           ' ----------------------------------'.format(VERSION)


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
        :param list_db:        list database entrys
        :type list_db:         boolean
        :param no_refine:      turn refinement off
        :type no_refine:       boolean
        :param invert:         invert the fragment during import, export or LS-fit
        :type invert:          boolean
        '''
        # options from the commandline options parser:
        self.options = OptionsParser(progname)
        self.external = False

        if not res_file_name and not external_restr:
            if not self.options.res_file:
                self.res_file = self.options.external_restr
                self.external = True
            else:
                self.res_file = self.options.res_file
                self.external = False
        else:
            if not res_file_name:
                self.res_file = external_restr
                self.external = True
            else:
                self.res_file = res_file_name
                self.external = False
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

        #  List of Database Fragments:
        if self.list_db:
            self.list_dbentrys()
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
        export = Export(self.export_fragment, self.invert)
        export.write_res_file()
        sys.exit(1)

    def export_to_clip(self):
        '''
        Exports the current fragment to the clipboard.
        '''
        from export import Export
        export = Export(self.export_clip)
        export.export_to_clip()
        sys.exit(1)

    def list_dbentrys(self):
        '''
        list all entries in the db.
        '''
        try:
            from terminalsize import get_terminal_size
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
        print('\n {} Fragments in the Database(s).'.format(num),
              '\n Feel free to add more fragments to "{}dsr_user_db.txt"' \
              '\n or mail them to dkratzert@gmx.de.'.format(dbdir + os.path.sep))

        for fragment in list(db.keys()):
            gdb.check_consistency(db[fragment], fragment)
            gdb.check_db_atom_consistency(db[fragment]['atoms'], fragment)
            gdb.check_db_header_consistency(db[fragment]['head'], fragment)
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
            names_list.append([i, fragname])
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
        print('\n\n Found following database entrys:\n')
        print(' Fragment          |  Full name, Comments')
        print(' ---------------------------------------------------------------------------')
        for line in results:
            print(' {:15s}   | {:20s}'.format(line[0], line[1]))
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
        shx.remove_afix()   # removes the afix 9


    def export_all_fragments(self):
        '''
        export all database entries at once
        '''
        from export import Export
        gdb = global_DB(self.invert)
        db = gdb.build_db_dict()
        dbnames = list(db.keys())
        for name in dbnames:
            export = Export(name)
            export.write_res_file()
        sys.exit(1)


    def set_final_db_sfac_types(self, db_atom_types, dbatoms, sfac_table):
        '''
        corrects the sfac types of the dbentry according to sfac card of the
        res file
        :param db_atom_types: element names of each atom in the database entry
                              like ['C', 'C', 'N', ... ]
        :type db_atom_types: list
        :param dbatoms: full atoms of the database entry
        :type dbatoms: list
        :param sfac_table: list of scattering factors from SHELXL
        :type sfac_table: list
        '''
        e2s = Elem_2_Sfac(sfac_table)
        atype = list(reversed(db_atom_types))
        for line in dbatoms:                    # go through db entry
            # replace scattering factor (line[1]) with true one
            line[1] = e2s.elem_2_sfac(atype.pop())
        return dbatoms


    def replacemode(self, res_target_atoms, rle, reslist, sfac_table):
        '''
        Target atoms are being replaced if this is executed
        '''
        for i in res_target_atoms:
            if '_' in i:
                print('\nDo you really want to REPLACE atom {} inside a residue?'.format(i))
                print('This will very likely not work.\n')
                break
        fa = FindAtoms(reslist)
        print('Replace mode active.')
        target_lines = fa.get_atom_line_numbers(res_target_atoms)
        #rp = False
        for i in target_lines:
            i = int(i)
            rle.remove_line(i, rem=False, remove=False, frontspace=True)
            #if reslist[i-1][:4] in ['RESI', 'PART']:
            #    rp = True
            #if reslist[i+1][:4] in ['RESI', 'PART'] and rp:
            #    print('empty resi/part?')
        fa.remove_adjacent_hydrogens(res_target_atoms, sfac_table)


    def go_refine(self, shx):
        '''
        actually starts the fragment fit
        '''
        try:
            shx.run_shelxl()
        except() as e:
            print(e)
            sys.exit()


    def generate_dfix_restraints(self, lf, reslist, dbatoms, residue_number, cell, part=''):
        '''
        returns a string of DFIX restraints for all 1,2- and 1,3-Bond distances
        in the current fragment.
        'DFIX at1 at2 distance\n DFIX at1 at2 distance\n ...'

        dbatoms: formated atoms (with new number scheme) of the fragment
        '''
        #from misc import format_atom_names
        #fa = FindAtoms(reslist)
        #lst_file = lf.read_lst_file()
        lst_file_coordinates = lf.get_all_coordinates
        lst_file_connectivity_table = lf.read_conntable()
        fragment_atoms = format_atom_names(dbatoms, part, residue_number)
        am = Adjacency_Matrix(fragment_atoms, lst_file_connectivity_table, lst_file_coordinates, cell)
        re = Restraints(lst_file_coordinates, am.get_adjmatrix, fragment_atoms, cell)
        dfixes = re.get_formated_12_dfixes+re.get_formated_13_dfixes+re.get_formated_flats
        return ''.join(dfixes)




    def main(self):
        '''
        main object to run DSR as command line program
        '''
        print(progname)
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

        fragment = dsrp.fragment

        fragline = gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
        dbatoms = gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
        dbhead = gdb.get_head_from_fragment(fragment)        # this is only executed once
        residue_class = gdb.get_resi_from_fragment(fragment)
        db_atom_types = get_atomtypes(dbatoms)                 # the atomtypes of the dbentry as list e.g. ['C', 'N', ...]
        resi = Resi(reslist, dsr_dict, dbhead, residue_class, find_atoms)
        dbhead = resi.make_resihead()
        sf = SfacTable(reslist, db_atom_types)
        sfac_table = sf.set_sfac_table()                 # from now on this sfac table is set

        ### corrects the atom type according to the previous defined global sfac table:
        dbatoms = self.set_final_db_sfac_types(db_atom_types, dbatoms, sfac_table)

        ## Insert FRAG ... FEND entry:
        rle.insert_frag_fend_entry(dbatoms, fragline, fvarlines)
        if not dsr_dict['resi'] and self.external:
            print('External restraint files are only possible in combination with residues!')
            sys.exit()

        print('Inserting {} into res File.'.format(fragment))
        if self.invert:
            print('Fragment inverted.')
        print('Source atoms:', ', '.join(dsrp.source))
        print('Target atoms:', ', '.join(dsrp.target))

        # several checks if the atoms in the dsr command line are consistent
        check_source_target(dsrp.source, dsrp.target, dbatoms)
        basefilename = filename_wo_ending(self.res_file)
        num = NumberScheme(reslist, dbatoms, resi.get_resinumber)
        fragment_numberscheme = num.get_fragment_number_scheme()
        afix = InsertAfix(reslist, dbatoms, db_atom_types, dbhead, dsr_dict,
                          sfac_table, find_atoms, fragment_numberscheme)
        afix_entry = afix.build_afix_entry(self.external, basefilename+'.dfix',
                                           resi.get_resiclass)
        # line where the dsr command is found in the resfile:
        dsr_line_number = dsrp.find_dsr_command(line=False)
        reslist[dsr_line_number] = reslist[dsr_line_number]+'\n'
        reslist.insert(dsr_line_number+1, afix_entry)

        ##### comment out all target atom lines in replace mode:
        if dsrp.command == 'REPLACE':
            self.replacemode(dsrp.target, rle, reslist, sfac_table)

        # write to file:
        shx = ShelxlRefine(reslist, basefilename, find_atoms)
        if not self.no_refine:
            shx.set_refinement_cycles('0')
        shx.remove_acta_card()

        rl.write_resfile(reslist, '.ins')

        #  Refine with L.S. 0 to insert the fragment
        if not self.no_refine:
            self.go_refine(shx)

        # Display the results from the list file:
        lf = ListFile(basefilename)
        lst_file = lf.read_lst_file()
        lfd = Lst_Deviations(lst_file)
        lfd.print_LS_fit_deviations()
        cell = lf.get_lst_cell_parameters

        # open res file again to restore 10 refinement cycles:
        rl = ResList(self.res_file)
        reslist = rl.get_res_list()
        shx = ShelxlRefine(reslist, basefilename, find_atoms)

        if not self.no_refine:
            self.set_post_refine_cycles(shx, '8')

        if dsrp.dfix_active:
            resinumber = resi.get_resinumber
            dfix_restraints = self.generate_dfix_restraints(lf,
                                            reslist,
                                            fragment_numberscheme,
                                            resinumber,
                                            cell,
                                            dsrp.part)
            if resinumber:
                if self.external:
                    external_file_name = write_dbhead_to_file(basefilename+'.dfix', dfix_restraints,
                                         resi.get_resiclass, resinumber)
                for n, line in enumerate(reslist):
                    if line.upper().startswith('RESI'):
                        if line.split()[1] == str(resinumber):
                            if self.external:
                                line = '{}\nREM The restraints for residue {} '\
                                'are in this file:\n+{}\n'.format(line, resinumber, external_file_name)
                            else:
                                line = '{} \n{}\n'.format(line, dfix_restraints)
                            reslist[n] = line
            else:
                line = '{}'.format(dfix_restraints) # insert restraints after dsr_line_number
                reslist[dsr_line_number-2] = reslist[dsr_line_number-2]+line

        if not self.no_refine:
            rl.write_resfile(reslist, '.res')



if __name__ == '__main__':
    '''main function'''
    #dsr = DSR(list_db=True)
    #dsr = DSR(res_file='p21c.res')
    dsr = DSR()



