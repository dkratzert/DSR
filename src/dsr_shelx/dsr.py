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

import os
import sys
from datetime import datetime
from pathlib import Path

import atomhandling
from afix import remove_duplicate_restraints, write_dbhead_to_file
from atomhandling import Elem_2_Sfac, rename_restraints_atoms
from constants import isoatomstr, sep_line, width
from dbfile import ImportGRADE, print_search_results, ParseDB, search_fragment_name
from dsrparse import DSRParser
from fit.quatfit import centroid, fit_fragment
from misc import cart_to_frac, chunks, frac_to_cart, touch, wrap_headlines, matrix_plus_vect, subtract_vect
from options import OptionsParser
from refine import ShelxlRefine
from resfile import ResList, ResListEdit, filename_wo_ending
from resi import Resi, remove_resi
from restraints import Restraints
from terminalsize import get_terminal_size
from version import VERSION
import time

minuse = ((width // 2) - 7) * '-'
program_name = '\n{} D S R - v{} {}'.format(minuse, VERSION, minuse)


class DSR(object):
    """
    main class
    """

    def __init__(self, options):
        """
        """
        time1 = time.perf_counter()
        # options from the commandline options parser:
        self.options = options
        self.external = False
        self.fragment = ''
        self.helpmsg = "*** Please ask dkratzert@gmx.de for help ***"
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
            import selfupdate
            selfupdate.update_dsr()
            sys.exit()
        #################################
        self.maindb_path = ''
        self.userdb_path = ''
        try:
            self.set_database_locations()
            self.gdb = ParseDB(self.maindb_path, self.userdb_path)
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
            raise
            sys.exit()
        #################################
        if not any(list(vars(self.options.all_options).values()) + [self.res_file]):
            self.options.error()
        if self.head_csv:
            self.head_to_gui()
        if self.search_extern:
            result = search_fragment_name(self.search_extern, self.gdb, numresults=8)
            for i in result:
                print(f'{i[0]};;{i[1]};;{i[2]};;{i[3]}')
            sys.exit()
        print(program_name)
        if self.list_db:
            self.list_dbentries()
        if self.search_string:
            result = search_fragment_name(self.search_string, self.gdb, numresults=7)
            print_search_results(result)
            sys.exit()
        # Export !all! fragments
        if self.export_all:
            self.export.export_all_fragments()
        # Export one fragment
        if self.export_fragment:
            print(f'Exporting "{self.fragment}" to {self.fragment}.res')
            self.fragment = self.export_fragment
            self.export.write_res_file(self.fragment)
            sys.exit()
        if self.export_clip:
            self.export.export_to_clip(self.fragment)
            sys.exit()
        # Import a GRADE fragment
        if self.import_grade:
            mog = ImportGRADE(self.import_grade, self.gdb, self.invert, self.maindb_path, self.userdb_path)
            mog.write_user_database()
            sys.exit()
        if not self.res_file:
            self.options.error()
        self.rl = ResList(self.res_file)
        self.reslist = self.rl.get_res_list()
        self.main()
        time2 = time.perf_counter()
        runtime = (time2 - time1)
        print('Runtime: {:>.1f} s'.format(runtime))
        print('DSR run complete.')

    ###############################################################################

    def set_database_locations(self):
        homedir = os.path.expanduser("~")
        self.userdb_path = os.path.join(homedir, "dsr_user_db.txt")
        if not os.path.isfile(self.userdb_path):
            touch(self.userdb_path)
        main_dbdir = os.environ.get("DSR_DIR", None)
        if main_dbdir:
            self.maindb_path = os.path.join(main_dbdir, 'dsr_db.txt')
        else:
            self.maindb_path = str(Path(__file__).parent / 'dsr_db.txt')
        return True

    def head_to_gui(self):
        """
        Exports current fragment header and atoms to the GUI
        """
        atoms = []
        try:
            atoms = self.export.export_to_gui(self.fragment)
        except Exception as e:
            # print(e)
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
        dbdir = os.path.expanduser('~')
        fragnames = []
        num = 0
        try:
            (width, _) = get_terminal_size()
        except():
            width = 80
        print('\n Entries found in the databases:\n')
        print(' Fragment         | Line | DB Name    | Full name, Comments ')
        print(sep_line)
        for num, line in enumerate(self.gdb.list_fragments()):
            fragnames.append(line[0])
            line = ' {:<17}| {:<5}| {:<11}| {}'.format(*line)
            print(line[:width - 1])
        print(f'\n {num} Fragments in the database(s).',
              f'\n Feel free to add more fragments to "{dbdir + os.path.sep}dsr_user_db.txt"\n')
        for fragment in fragnames:
            self.gdb.check_consistency(fragment)
            self.gdb.check_db_atom_consistency(fragment)
            self.gdb.check_db_restraints_consistency(fragment)
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
        basefilename = filename_wo_ending(self.res_file)
        if not basefilename:
            print('*** Illegal option ***')
            sys.exit()
        if len(self.reslist) == 0:
            print("*** The input file is empty. Can not proceed! ***")
            sys.exit()
        find_atoms = atomhandling.FindAtoms(self.reslist)
        rle = ResListEdit(self.reslist, find_atoms)
        dsrp = DSRParser(self.reslist)
        self.fragment = dsrp.fragment
        restraints = self.gdb.get_restraints(self.fragment)  # this is only executed once
        db_residue_string = self.gdb.get_resi(self.fragment)
        dbatoms = self.gdb.get_atoms(self.fragment, self.invert)  # only the atoms of the dbentry as list
        # the atomtypes of the dbentry as list e.g. ['C', 'N', ...]
        db_atom_types = atomhandling.get_atomtypes(dbatoms)
        sf = atomhandling.SfacTable(self.reslist, db_atom_types)
        sfac_table = sf.set_sfac_table()  # from now on this sfac table is set
        resi = Resi(dsrp, db_residue_string, find_atoms)
        # line where the dsr command is found in the resfile:
        if dsrp.cf3_active:
            from cf3fit import CF3
            cf3 = CF3(rle, find_atoms, self.reslist, self.fragment, sfac_table,
                      basefilename, dsrp, resi, self.res_file, self.options)
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
        self.gdb.check_db_restraints_consistency(self.fragment)
        self.gdb.check_sadi_consistence(self.fragment)
        if dsrp.occupancy:
            rle.set_free_variables(dsrp.occupancy)
        restraints = remove_resi(restraints)
        # corrects the atom type according to the previous defined global sfac table:
        dbatoms = atomhandling.set_final_db_sfac_types(db_atom_types, dbatoms, sfac_table)
        if not dsrp.unit_line:
            print('*** No UNIT instruction in res file found! Can not proceed! ***')
        print('Inserting {} into res File.'.format(self.fragment))
        if self.invert:
            print('Fragment inverted.')
        print('Source atoms: {}'.format(', '.join(dsrp.source)))
        print('Target atoms: {}'.format(', '.join(dsrp.target)))
        shx = ShelxlRefine(self.reslist, basefilename, find_atoms, self.options)
        shx.backup_shx_file()
        # several checks if the atoms in the dsr command line are consistent
        atomhandling.check_source_target(dsrp.source, dsrp.target, dbatoms)
        num = atomhandling.NumberScheme(self.reslist, dbatoms, dsrp)
        # returns also the atom names if residue is active
        fragment_numberscheme = num.get_fragment_number_scheme()
        print('Fragment atom names: {}'.format(', '.join(fragment_numberscheme)))
        dfix_head = ''
        if dsrp.dfix:
            restr = Restraints(self.fragment, self.gdb)
            dfix_12 = restr.get_formated_12_dfixes()
            dfix_13 = restr.get_formated_13_dfixes()
            flats = restr.get_formated_flats()
            restraints = dfix_12 + dfix_13 + flats
        # ##########Not using SHELXL for fragment fit: ###########

        print("--- Using fast fragment fit ---")
        if self.options.target_coords:
            target_coords = chunks(self.options.target_coords, 3)
        else:
            # {'C1': ['1.123', '0.7456', '3.245']}
            target_coordinates = find_atoms.get_atomcoordinates(dsrp.target)
            target_coords = [target_coordinates[key] for key in dsrp.target]
        # Uppercase is important here to avoid KeyErrors in source_atoms generation
        atnames = self.gdb.get_atomnames(self.fragment, uppercase=True)
        source_atoms = dict(zip(atnames, self.gdb.get_coordinates(self.fragment, cartesian=True, invert=self.invert)))
        # Coordinates only from the source, not the entire fragment:
        source_coords = [source_atoms[x] for x in dsrp.source]
        target_coords = [frac_to_cart(x, rle.cell) for x in target_coords]
        # The source and target atom coordinates are fitted first. Then The complete fragment
        # is rotated and translated to the target position as calculated before.
        # parameter cartiesian has to be false here:
        fragment_coords = self.gdb.get_coordinates(self.fragment, cartesian=True, invert=self.invert)
        fitted_fragment, rmsd = fit_fragment(fragment_coords,
                                             source_atoms=source_coords,
                                             target_atoms=target_coords)
        # Moving back to the position of the first atom to have a reference:
        # I have to make sure that I use the centroid of the correct atoms from target and source,
        # otherwise the fragment is shifted to a wrong position.
        # The third atom from the fragment e.g. has to be the third from the fragment to get
        # the correct centroid:
        center_difference = subtract_vect(centroid(target_coords),
                                          centroid([list(fitted_fragment)[atnames.index(dsrp.source[x])] for x in
                                                    range(len(source_coords))]))
        # finishing shift to correct centroid:
        fitted_fragment = matrix_plus_vect(fitted_fragment, center_difference)
        # Or even lower than 0.1?
        if rmsd < 0.1:
            print('Fragment fit successful with RMSD of: {:8.3}'.format(rmsd))
        else:
            print('*** Fragment fit might have failed with RMSD of: {:8.3} ***'.format(rmsd))
        fitted_fragment = [cart_to_frac(x, rle.cell) for x in fitted_fragment]
        afix_entry = []
        e2s = Elem_2_Sfac(sfac_table)
        for at, coord, atype in zip(fragment_numberscheme, fitted_fragment, db_atom_types):
            sfac_num = str(e2s.elem_2_sfac(atype))
            if dsrp.occupancy:
                occ = float(dsrp.occupancy)
            else:
                occ = 11.0
            afix_entry.append(isoatomstr.format(at, sfac_num, coord[0], coord[1], coord[2], occ, 0.03))
        afix_entry = "\n".join(afix_entry)
        new_atomnames = list(reversed(fragment_numberscheme))
        same_resi = ''
        if not dsrp.resiflag:
            restraints = rename_restraints_atoms(new_atomnames, self.gdb.get_atomnames(self.fragment), restraints)
        else:
            restraints = resi.format_restraints(restraints)
            # SADI\n
            same_resi = ["SAME_{} {} > {}\n".format(resi.get_residue_class, new_atomnames[-1], new_atomnames[0])]
        # Adds a "SAME_resiclass firstatom > lastatom" to the afix:
        if not self.options.rigid_group:
            restraints += same_resi
            # if dsrp.resiflag:  # <- Or should I do this?
            restraints += ["SIMU 0.04 0.08 1"]
        if not self.options.external_restr:
            restraints = remove_duplicate_restraints(self.reslist, restraints, resi.get_residue_class)
        restraints = wrap_headlines(restraints)
        dfx_file_name = ''
        if dsrp.part:
            afix_entry = "PART {}  {}\n".format(dsrp.part, dsrp.occupancy) + afix_entry + "\nPART 0"
        if dsrp.resiflag:
            afix_entry = 'RESI {} {}\n{}\nRESI 0'.format(resi.get_residue_class, resi.get_resinumber, afix_entry)
        if self.options.rigid_group:
            afix_entry = 'AFIX 9\n' + afix_entry
        if self.options.external_restr and not self.rigid:
            pname, ext = os.path.splitext(basefilename + '.dfix')
            if dsrp.dfix:
                dfx_file_name = pname + "_dfx" + ext
            else:
                dfx_file_name = pname + ext
            dfx_file_name = write_dbhead_to_file(dsrp, dfx_file_name, restraints, resi.get_residue_class,
                                                 resi.get_resinumber)
            if dsrp.resiflag:
                restraints = 'REM Restraints for residue {}:\n+{}\n'.format(resi.get_residue_class, dfx_file_name)
            else:
                restraints = 'REM Restraints for DSR fragment:\n+{}\n'.format(dfx_file_name)
        if self.options.rigid_group:
            afix_entry += '\nAFIX 0\n'

        # Adds the origin of restraints and fragment to res file:
        import textwrap
        source = textwrap.wrap("REM Restraints for Fragment {}, {} from: {}. "
                               "Please cite https://doi.org/10.1107/S1600576718004508".format(
            self.fragment,
            self.gdb.get_fragment_name(self.fragment),
            self.gdb.get_src(self.fragment)),
            width=74, subsequent_indent='REM ')
        if dsrp.resi:
            hfixes = '\n'.join(self.gdb.get_hfixes(self.fragment, resi.get_residue_class))
            if hfixes:
                hfixes += '\n'
        else:
            hfixes = ''
        source = '\n'.join(source) + '\n' + hfixes
        # check if restraints already inserted:
        for line in self.reslist:
            try:
                if line.split()[4] == self.fragment + ',':
                    source = ''
                    break
            except IndexError:
                continue
        # + 'AFIX 0\n' before hklf seems to be not needed after shelx-2013:
        self.reslist[dsrp.hklf_line - 1] = self.reslist[dsrp.hklf_line - 1] + afix_entry + '\n'
        if not self.rigid:
            self.reslist[dsrp.unit_line] = self.reslist[dsrp.unit_line] + source + ''.join(restraints)

        # write to file:
        self.rl.write_resfile(self.reslist, '.res')
        if dsrp.command == 'REPLACE':
            print("Replace mode active\n")
            self.rl = ResList(self.res_file)
            reslist = self.rl.get_res_list()
            self.reslist, find_atoms = atomhandling.replace_after_fit(self.rl, reslist, resi,
                                                                      fragment_numberscheme, rle.get_cell())
            self.rl.write_resfile(self.reslist, '.res')
        os.remove(shx.backup_file)


class Multilog(object):
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
            res = ''
            for f in self._files:
                res = getattr(f, attr, *args)(*a, **kw)
            return res

        return g


def main():
    options = None
    lstfile = ''
    is_listfile = False
    try:
        lstfile = open('./dsr-log.lst', 'w')
    except IOError:
        pass
    else:
        sys.stdout = Multilog([sys.stdout, lstfile])
        sys.stderr = Multilog([sys.stderr, lstfile])
        is_listfile = True
    try:
        options = OptionsParser(program_name)
        if is_listfile:
            lstfile.write('Python version: {}\n'.format(sys.version))
            lstfile.write("Date: {} \n".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            lstfile.write(options.__str__())
        DSR(options)
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
        if options:
            print('Commandline: {}'.format(options.all_options))
        print('Platform: {} {}, {}'.format(platform.system(),
                                           platform.release(), ' '.join(platform.uname())))
        print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        raise


if __name__ == '__main__':
    main()
