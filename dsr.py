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
from dbfile import global_DB, search_fragment_name
from constants import width, sep_line
from misc import reportlog, remove_file, find_line, remove_line
from options import OptionsParser
from os.path import expanduser
from terminalsize import get_terminal_size
from dsrparse import DSR_Parser
from dbfile import ImportGRADE, print_search_results
from resi import Resi
from restraints import ListFile, Lst_Deviations, Restraints
from afix import Afix
from refine import ShelxlRefine
from resfile import ResList, filename_wo_ending, ResListEdit

VERSION = '203'
# dont forget to change version in Innoscript file, spec file and deb file.

program_name = '\n'+((width//2)-9)*'-'+' D S R - v{} '.format(VERSION)+((width//2)-8)*'-'

# TODO and ideas:
"""
- Add Rcomplete

From SHELXL user guide:
A free variable is a refinable parameter that can be used to impose a variety of additional
linear constraints, e.g. to atomic coordinates, occupancies or displacement parameters.
Starting values for all free variables are supplied on the FVAR instruction. Since the first
FVAR parameter is the overall scale factor, there is no free variable number 1. If an atom
parameter is given a value greater than 15 or less than -15, it is interpreted as a reference to
a free variable. A positive value (10k+p) is decoded as p times free variable number k [fv(k)],
and a negative value (i.e. k and p both negative) means p times [fv(–k)–1].

- port to JANA?
  -> learn JANA
  -> What do I need to change?

"""


class DSR():
    """
    main class
    """
    def __init__(self, options):
        """
        """
        import time
        time1 = time.clock()
        # options from the commandline options parser:
        self.options = options
        self.external = False
        self.fragment = ''
        self.helpmsg = "*** Please ask daniel.kratzert@ac.uni-freiburg.de for help ***"
        self.res_file = self.options.res_file
        if self.options.external_restr:
            self.external = True
            self.res_file = self.options.external_restr
        self.export_fragment = self.options.export_fragment
        if self.export_fragment:
            self.res_file = self.options.export_fragment
        if self.export_fragment:
            self.fragment = self.export_fragment
        self.export_clip = self.options.export_clip
        if self.export_clip:
            self.fragment = self.export_clip
        self.import_grade = self.options.import_grade
        self.export_all = self.options.export_all
        self.list_db = self.options.list_db
        self.list_db_csv = self.options.list_db_csv
        self.no_refine = self.options.no_refine
        self.invert = self.options.invert
        self.rigid = self.options.rigid_group
        self.search_string = self.options.search_string
        self.search_extern = self.options.search_extern
        self.head_csv = self.options.head_for_gui
        if self.head_csv:
            self.fragment = self.head_csv
        if self.options.selfupdate:
            print(program_name)
            import selfupdate
            selfupdate.update_dsr()
            sys.exit()
        #################################
        try:
            self.gdb = global_DB(invert=self.invert)
        except Exception as e:  # @UnusedVariable
            print("*** Initializing the database failed ***")
            raise
        #  List of database Fragments:
        if self.list_db_csv:
            print('DSR version: {}'.format(VERSION))
            for i in self.gdb.list_fragments():
                print('{};;{};;{};;{}'.format(i[0], i[3], i[1], i[2]))
            sys.exit()
        try:
            from export import Export
            self.export = Export(gdb=self.gdb, invert=self.invert)
        except Exception as e:
            print("*** Unable to export informations from DSR ***")
            print(e)
            sys.exit()
        #################################
        if self.head_csv:
            self.head_to_gui()
        if self.search_extern:
            result = search_fragment_name(self.search_extern, self.gdb, numresults=7)
            for i in result:
                print('{};;{};;{};;{}'.format(i[0], i[1], i[2], i[3]))
            sys.exit()
        print(program_name)
        if self.list_db:
            self.list_dbentries()
        if self.search_string:
            result = search_fragment_name(self.search_string, self.gdb)
            print_search_results(result)
            sys.exit()
        # Export !all! fragments
        if self.export_all:
            self.export.export_all_fragments()
        # Export one fragment
        if self.export_fragment:
            print('Exporting "{0}" to {0}.res'.format(self.fragment))
            self.fragment = self.export_fragment
            try:
                self.export.write_res_file(self.fragment)
            except:
                raise
            sys.exit()
        if self.export_clip:
            try:
                self.export.export_to_clip(self.fragment)
            except:
                raise
            sys.exit()
        # Import a GRADE fragment
        if self.import_grade:
            mog = ImportGRADE(self.import_grade, self.invert)
            mog.write_user_database()
            sys.exit()
        if not any(list(vars(self.options.all_options).values())+[self.res_file]):
            self.options.error()
        if self.res_file == False:
            self.options.error()
        self.rl = ResList(self.res_file)
        self.reslist = self.rl.get_res_list()
        self.main()
        time2 = time.clock()
        runtime = (time2 - time1)
        print('Runtime: {:>.1f} s'.format(runtime))
        print('DSR run complete.')

###############################################################################

    def head_to_gui(self):
        """
        Exports current fragment header and atoms to the GUI
        """
        atoms = []
        try:
            atoms = self.export.export_to_gui(self.fragment)
        except:
            print("*** Could not get atom information ***")
            print(self.helpmsg)
            sys.exit()
        print("\n<atoms>")
        print(atoms)
        print("</atoms>")
        # prints most of the needed info:
        self.gdb.get_head_for_gui(self.fragment)
        sys.exit()

    def list_dbentries(self):
        """
        list all entries in the db.
        """
        dbdir = expanduser('~')
        fragnames = []
        num = 0
        try:
            (width, height) = get_terminal_size()  # @UnusedVariable
        except():
            width = 80
        print('\n Entries found in the databases:\n')
        print(' Fragment         | Line | DB Name    | Full name, Comments ')
        print(sep_line)
        for num, line in enumerate(self.gdb.list_fragments()):
            fragnames.append(line[0])
            line = ' {:<17}| {:<5}| {:<11}| {}'.format(*line)
            print(line[:width - 1])
        print('\n {} Fragments in the database(s).'.format(num),
              '\n Feel free to add more fragments to "{}dsr_user_db.txt"'.format(dbdir + os.path.sep))
        for fragment in fragnames:
            self.gdb.check_consistency(fragment)
            self.gdb.check_db_atom_consistency(fragment)
            self.gdb.check_db_header_consistency(fragment)
            self.gdb.check_sadi_consistence(fragment)
        from selfupdate import is_update_needed
        if is_update_needed(silent=True):
            print("\n*** An update for DSR is available. You can update with 'dsr -u' ***")
        sys.exit()

    def main(self):
        """
        main object to run DSR as command line program
        """
        dbatoms = []
        # The database content:
        import atomhandling
        basefilename = filename_wo_ending(self.res_file)
        if not basefilename:
            print('*** Illegal option ***')
            sys.exit()
        if len(self.reslist) == 0:
            print("*** The input file is empty. Can not proceed! ***")
            sys.exit()
        find_atoms = atomhandling.FindAtoms(self.reslist)
        rle = ResListEdit(self.reslist, find_atoms)
        dsrp = DSR_Parser(self.reslist, rle)
        dsr_dict = dsrp.get_dsr_dict
        fvarlines = rle.find_fvarlines()
        self.fragment = dsrp.fragment.lower()
        dbhead = self.gdb.get_head_from_fragment(self.fragment)        # this is only executed once
        db_residue_string = self.gdb.get_resi_from_fragment(self.fragment)
        dbatoms = self.gdb.get_atoms_from_fragment(self.fragment)      # only the atoms of the dbentry as list
        # the atomtypes of the dbentry as list e.g. ['C', 'N', ...]
        db_atom_types = atomhandling.get_atomtypes(dbatoms)
        sf = atomhandling.SfacTable(self.reslist, db_atom_types)
        sfac_table = sf.set_sfac_table()                 # from now on this sfac table is set
        resi = Resi(self.reslist, dsr_dict, dbhead, db_residue_string, find_atoms)
        # line where the dsr command is found in the resfile:
        dsr_line_number = dsrp.find_dsr_command(line=False)
        if dsrp.cf3_active:
            from cf3fit import CF3
            cf3 = CF3(rle, find_atoms, self.reslist, self.fragment, sfac_table,
                      basefilename, dsr_dict, resi, self.res_file, self.options)
            if self.fragment == 'cf3':
                cf3.cf3(afix='130')
            if self.fragment == 'cf6':
                cf3.cf3(afix='120')
            if self.fragment == 'cf9':
                cf3.cf9()
            print('\nFinished...')
            sys.exit()
        # checks have to be after CF3, CF6 etc.
        self.gdb.check_consistency(self.fragment)
        self.gdb.check_db_atom_consistency(self.fragment)
        self.gdb.check_db_header_consistency(self.fragment)
        self.gdb.check_sadi_consistence(self.fragment)
        if dsrp.occupancy:
            rle.set_free_variables(dsrp.occupancy)
        fragline = self.gdb.get_fragline_from_fragment(self.fragment)  # full string of FRAG line
        dbhead = resi.remove_resi(dbhead)
        # corrects the atom type according to the previous defined global sfac table:
        dbatoms = atomhandling.set_final_db_sfac_types(db_atom_types, dbatoms, sfac_table)

        # Insert FRAG ... FEND entry:
        rle.insert_frag_fend_entry(dbatoms, fragline, fvarlines)

        print('Inserting {} into res File.'.format(self.fragment))
        if self.invert:
            print('Fragment inverted.')
        print('Source atoms: {}'.format(', '.join(dsrp.source)))
        print('Target atoms: {}'.format(', '.join(dsrp.target)))

        # several checks if the atoms in the dsr command line are consistent
        atomhandling.check_source_target(dsrp.source, dsrp.target, dbatoms)
        num = atomhandling.NumberScheme(self.reslist, dbatoms, resi.get_resinumber)
        # returns also the atom names if residue is active
        fragment_numberscheme = num.get_fragment_number_scheme()
        print('Fragment atom names: {}'.format(', '.join(fragment_numberscheme)))
        dfix_head = ''
        if dsrp.dfix_active:
            restr = Restraints(self.export, self.fragment, self.gdb)
            dfix_12 = restr.get_formated_12_dfixes()
            dfix_13 = restr.get_formated_13_dfixes()
            flats = restr.get_formated_flats()
            dfix_head = dfix_12+dfix_13+flats
        afix = Afix(self.reslist, dbatoms, db_atom_types, dbhead, dsr_dict,
                    sfac_table, find_atoms, fragment_numberscheme, self.options, dfix_head)
        afix_entry = afix.build_afix_entry(self.external, basefilename+'.dfix', resi)
        if dsr_line_number < fvarlines[-1]:
            print('\n*** Warning! The DSR command line MUST NOT appear before FVAR '
                  'or the first atom in the .res file! ***')
            print('*** Can not proceed... ***\n')
            sys.exit()
        # Adds the origin of restraints and fragment to res file:
        import textwrap
        source = textwrap.wrap("REM Restraints for Fragment {}, {} from: {}. "
                               "Please cite doi:10.1107/S1600576715005580".format(
                                    self.fragment,
                                    self.gdb.get_name_from_fragment(self.fragment),
                                    self.gdb.get_src_from_fragment(self.fragment)),
                                width=74, subsequent_indent='REM ')
        # TODO: test if slow for big files:
        for line in self.reslist:
            try:
                if line.split()[4] == self.fragment + ',':
                    source = ''
                    break
            except IndexError:
                continue
        self.reslist[dsr_line_number] = self.reslist[dsr_line_number] + '\n' + '\n'.join(source) + '\n'+ afix_entry
        # write to file:
        shx = ShelxlRefine(self.reslist, basefilename, find_atoms, self.options)
        acta_lines = shx.remove_acta_card()
        cycles = shx.get_refinement_cycles
        shx.set_refinement_cycles('0')  # also sets shx.cycles to current value
        self.rl.write_resfile(self.reslist, '.ins')
        if self.no_refine:
            print('\nPerforming no fragment fit. Just prepared the .ins file for you.')
            return
        #  Refine with L.S. 0 to insert the fragment
        try:
            shx.run_shelxl()
        except:
            raise
            sys.exit()
        # Display the results from the list file:
        lf = ListFile(basefilename)
        lst_file = lf.read_lst_file()
        shx.check_refinement_results(lst_file)
        lfd = Lst_Deviations(lst_file)
        lfd.print_LS_fit_deviations()
        cell = rle.get_cell()
        # open res file again to restore 8 refinement cycles:
        self.rl = ResList(self.res_file)
        reslist = self.rl.get_res_list()
        # remove the "REM " instriction bevore the +dfixfile instruction
        plusline = find_line(reslist, "REM "+afix.rand_id_dfx)
        if plusline:
            reslist[plusline-1] = reslist[plusline-1][4:]
            remove_line(reslist, plusline, remove=True)
        if dsrp.command == 'REPLACE':
            reslist, find_atoms = atomhandling.replace_after_fit(self.rl, reslist, resi,
                                                    fragment_numberscheme, cell)
        shx = ShelxlRefine(reslist, basefilename, find_atoms, self.options)
        shx.restore_acta_card(acta_lines)
        try:
            if cycles:
                shx.set_refinement_cycles(cycles) # restores last LS value
            else:
                shx.set_refinement_cycles(8)
        except IndexError:
            print('*** Unable to set refinement cycles ***')
        if not self.options.rigid_group:
            shx.remove_afix(afix.rand_id_afix)   # removes the afix 9
        # final resfile write:
        self.rl.write_resfile(reslist, '.res')


class multilog(object):
    """
    This class copies all output from stdout and stderr to a file
    It acts like tee with following usage:
    sys.stdout = multifile([sys.stdout, lstfileobj])
    """
    def __init__(self, files):
        self._files = files

    def __getattr__(self, attr, *args):
        return self._wrap(attr, *args)

    def _wrap(self, attr, *args):
        def g(*a, **kw):
            for f in self._files:
                res = getattr(f, attr, *args)(*a, **kw)
            return res
        return g


if __name__ == '__main__':
    '''main function'''
    lstfile = ''
    is_listfile = False
    try:
        lstfile = open('./dsr-log.lst', 'w')
    except IOError:
        pass
    else:
        sys.stdout = multilog([sys.stdout, lstfile])
        sys.stderr = multilog([sys.stderr, lstfile])
        is_listfile = True
    try:
        options = OptionsParser(program_name)
        if is_listfile:
            lstfile.write('Python version: {}\n'.format(sys.version))
        dsr = DSR(options)
    except Exception:
        import platform
        if is_listfile:
            lstpath = os.path.abspath(lstfile.name)
            lst = 'the file "{}" \nand '.format(lstpath)
        else:
            lst = "this error message and "
        print('\n*** Congratulations! You found a bug in DSR. Please send {}the .res file '
              '(if possible) to dkratzert@gmx.de ***\n\n'.format(lst))
        print('DSR version: {}'.format(VERSION))
        print('Python version: {}'.format(sys.version))
        print('Platform: {} {}, {}'.format(platform.system(),
                                            platform.release(), ' '.join(platform.uname())))
        raise