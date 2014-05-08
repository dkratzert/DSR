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
from atomhandling import SfacTable, get_atomtypes, check_source_target, Elem_2_Sfac
from atomhandling import FindAtoms, NumberScheme
from afix import InsertAfix, write_dbhead_to_file
from refine import ShelxlRefine
from resi import Resi
from misc import get_replace_mode, find_line
from restraints import ListFile, Lst_Deviations, format_atom_names
from restraints import Restraints, Adjacency_Matrix

# TODO and ideas:
# cell during export too wide
# -automatically invert fragment if fit fails and fit again?
# -FLAT: see for every ring atom if there is also an atteched atom in the same plane
# -To the manual: Avogradro/res-file -> rename -> mercury -> mol2-file -> GRADE
# -To manual: what to do if something does not work.
# -try to iterate all target atoms if fit is really bad
# -select atoms in shelxle > export to db
# -If second part of a disordered atom is very close to the first, make EADP?
# -detect empty residues and parts after atom deletion
# -debian package: /usr/src/packages/BUILD # dpkg-deb --build dsr


VERSION = '1.4.6'
progname = '\n-----------------------------'\
           ' D S R - v{}' \
           ' ----------------------------------'.format(VERSION)

class DSR():
    '''
    main class 
    '''
    def __init__(self, res_file=None, external_restr=None, export_fragment=None, 
        export_clip=None, import_grade=None, export_all=None, list_db=None, no_refine=None, invert=None):
        # options from the commandline options parser:
        self.options = OptionsParser(progname)
        self.external = False
        
        if not res_file and not external_restr:
            if not self.options.res_file:
                self.res_file = self.options.external_restr
                self.external = True
            else:
                self.res_file = self.options.res_file
                self.external = False
        else:
            if not res_file:
                self.res_file = external_restr
                self.external = True
            else:
                self.res_file = res_file
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

        #  List of Database Fragments:   
        if self.list_db:
            self.list_dbentrys()
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
        runtime = (time2-time1)
        print('Runtime: {:>.1f} s'.format(runtime))
        print('\nDSR run complete.')
        

    def do_export_fragment(self):
        ''' 
        Exports the current fragment.
        '''
        from export import Export
        export = Export(self.export_fragment, self.invert)
        export.write_file()
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
            (width, height) = get_terminal_size()
        except():
            width = 80
        gdb = global_DB()
        db = gdb.build_db_dict()
        print('\n Entries found in the databases:\n')
        print(' Fragment         | Line | DB Name    | Full name, Comments ')
        print(' ---------------------------------------------------------------------------')
    
        frags = sorted(db.keys())
        for i in frags:
            line = ' {:<17}| {:<5}| {:<11}| {}'.format(
                    i, gdb.get_line_number_from_fragment(i), gdb.get_db_from_fragment(i), 
                    ', '.join(gdb.get_comment_from_fragment(i)))
            print(line[:width-1])
        try:
            if os.environ["DSR_DB_DIR"]:
                dbdir = os.environ["DSR_DB_DIR"]
        except(KeyError):
            dbdir = '.'
        print('\n Feel free to add more fragments to "{}dsr_user_db.txt"' \
              '\n or mail them to dkratzert@gmx.de.'.format(dbdir+os.path.sep))
        
        for fragment in list(db.keys()):
            gdb.check_consistency(db[fragment], fragment)
            gdb.check_db_atom_consistency(db[fragment]['atoms'], fragment)
            gdb.check_db_header_consistency(db[fragment]['head'], fragment)
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
            export.write_file()
        sys.exit(1)
    
    
    def set_final_db_sfac_types(self, dbtypes, dbatoms, sfac_table):
        '''
        corrects the sfac types of the dbentry according to sfac card of the
        res file
        '''
        e2s = Elem_2_Sfac(sfac_table)
        atype = list(reversed(dbtypes))
        for line in dbatoms:                             # go through db entry
            line[1] = e2s.elem_2_sfac(atype.pop())       # replace scattering factor (line[1]) with true one
    
    
    def replacemode(self, res_target_atoms, rle, reslist):
        '''
        Target atoms are being replaced if this is executed
        '''
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
        fa.remove_adjacent_hydrogens(res_target_atoms)    
    
    
    def go_refine(self, shx):
        '''
        actually starts the fragment fit
        '''
        try:
            shx.run_shelxl()
        except() as e:
            print(e)
            sys.exit() 

    
    def generate_dfix_restraints(self, lf, reslist, dbatoms, residue, cell, part=''):
        '''
        returns a string of DFIX restraints for all 1,2- and 1,3-Bond distances
        in the current fragment.
        'DFIX at1 at2 distance\n DFIX at1 at2 distance\n ...'
        
        dbatoms: formated atoms (with new number scheme) of the fragment
        '''
        from misc import format_atom_names
        fa = FindAtoms(reslist)
        lst_file = lf.read_lst_file()
        coords = lf.get_all_coordinates
        conntable = lf.read_conntable()
        fragment_atoms = format_atom_names(dbatoms, part, residue)
        am = Adjacency_Matrix(fragment_atoms, conntable, coords, cell)
        re = Restraints(coords, am.get_adjmatrix, fragment_atoms, cell)
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
            rle.set_fvar(dsrp.occupancy, fvarlines)
        
        fragment = dsrp.fragment
    
        fragline = gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
        dbatoms = gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
        dbhead = gdb.get_head_from_fragment(fragment)        # this is only executed once
        residue = gdb.get_resi_from_fragment(fragment)
        dbtypes = get_atomtypes(dbatoms)                 # the atomtypes of the dbentry as list e.g. ['C', 'N', ...]
        resi = Resi(reslist, dsr_dict, dbhead, residue, find_atoms)
        dbhead = resi.make_resihead()
        sf = SfacTable(reslist, dbtypes, self.res_file)
        sfac_table = sf.set_sfac_table()                 # from now on this sfac table is set
        
        ### corrects the atom type according to the previous defined global sfac table:
        self.set_final_db_sfac_types(dbtypes, dbatoms, sfac_table)
    
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
        numberscheme = num.get_fragment_number_scheme()
        afix = InsertAfix(reslist, dbatoms, dbtypes, dbhead, dsr_dict, 
                          sfac_table, find_atoms, numberscheme)
        afix_entry = afix.build_afix_entry(self.external, basefilename+'.dfx', resi.get_resiclass)
        # line where the dsr command is found in the resfile:
        dsrline = dsrp.find_dsr_command(line=False) 
        reslist[dsrline] = reslist[dsrline]+'\n'
        reslist.insert(dsrline+1, afix_entry)
    
        ##### comment out all target atom lines in replace mode:  
        if dsrp.command == 'REPLACE':
            self.replacemode(dsrp.target, rle, reslist)

        # write to file:
        shx = ShelxlRefine(reslist, basefilename, find_atoms)
        if not self.no_refine:
            shx.set_refinement_cycles('0')
        shx.remove_acta()
        
        rl.write_resfile(reslist, '.ins')  
        
        #  Refine with L.S. 0 to insert the fragment
        if not self.no_refine:
            self.go_refine(shx)
        
        # Display the results from the list file:
        lf = ListFile(basefilename)
        lst_file = lf.read_lst_file()
        lfd = Lst_Deviations(lst_file)
        lfd.print_deviations()
        cell = lf.get_cell_params
        
        # open res file again to restore 10 refinement cycles:
        rl = ResList(self.res_file)
        reslist = rl.get_res_list()
        shx = ShelxlRefine(reslist, basefilename, find_atoms)
    
        if not self.no_refine:
            self.set_post_refine_cycles(shx, '8')

        if dsrp.dfix:
            resinumber = resi.get_resinumber
            dfix = self.generate_dfix_restraints(lf, 
                                            reslist, 
                                            numberscheme, 
                                            resinumber,
                                            cell,
                                            dsrp.part)
            if resinumber:
                if self.external:
                    externalfile_name = write_dbhead_to_file(basefilename+'.dfx', dfix, resi.get_resiclass, resinumber)
                for n, line in enumerate(reslist):
                    if line.upper().startswith('RESI'):
                        if line.split()[1] == str(resinumber):
                            if self.external:
                                line = '{}\nREM The restraints for residue {} are in this file:\n+{}\n'.format(line, resinumber, externalfile_name) 
                            else:
                                line = '{} \n{}\n'.format(line, dfix) 
                            reslist[n] = line    
            else:
                line = '{}'.format(dfix) # insert restraints after dsrline
                reslist[dsrline-2] = reslist[dsrline-2]+line
 
        if not self.no_refine:
            rl.write_resfile(reslist, '.res')
        
  
    
if __name__ == '__main__':
    '''main function'''
    #dsr = DSR(no_refine=False, res_file='p21c.res')
    dsr = DSR()

    

