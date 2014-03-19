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
from afix import InsertAfix
from refine import ShelxlRefine
from resi import Resi
from misc import get_replace_mode, find_line
from restraints import ListFile, Lst_Deviations, format_atom_names
from restraints import Restraints, Adjacency_Matrix

# TODO and ideas:
# -check for residues with same class and differing atom names inside.
# -count number of atoms in old and new residues to compare them for error checking.
# -To the manual: Avogradro/res-file -> rename -> mercury -> mol2-file -> GRADE
# -To manual: what to do if something does not work.
# -more detailed comments for grade imports:
#   REM Produced by Grade Web Server http://grade.globalphasing.org
#   REM GEN: Generated by GRADE 1.2.5 (December 20 2013)
#   REM GEN: from SMILES C1=CC=C2C(=C1)C=C(C3=CC=CC=C23)Br
#   REM GEN: using MOGUL 1.6(RC5), CSD as535be with quantum mechanics RM1
#   REM grade-cif2shelx output
#   REM Version: 0.0.5 <Dec 20 2013>
# -try to iterate all target atoms if fit is really bad
# -Iteratively analyze disagreeable restraints and set ESDs of relative
#  restraints accordingly. (Only if shift is not too big, has to be nearly converged)
# -select atoms in shelxle > export to db
# -if part > 2 and same residue and positive occ then SUMP
# -use +filename.dfix for restraints.
# -If second part of a disordered atom is very close to the first, make EADP?
# -If e.g. OCC in dsrline and index(OCC)+1 != ' ', then dsrline.insert(' ', index(OCC)+1)
# -should I make rem dsr IMPORT from Atom1 to Atom2 ?
#  with or without their restraints?
# -detect empty residues and parts after atom deletion
# -add SIMU and RIGU after Grade import
# -debian package: /usr/src/packages/BUILD # dpkg-deb --build dsr
# -check atoms bond valency after fit to decide if fit was sucessful.

VERSION = '1.3.3'
progname = '\n-----------------------------'\
           ' D S R - v{}' \
           ' ----------------------------------'.format(VERSION)

class DSR():
    '''
    main class
    '''
    def __init__(self, res_file=None, export_fragment=None, import_grade=None, 
                export_all=None, list_db=None, no_refine=None):
        # options from the commandline options parser:
        self.options = OptionsParser(progname)
        if not res_file:
            self.res_file = self.options.res_file
        else:
            self.res_file = res_file
        if not export_fragment:
            self.export_fragment = self.options.export_fragment
        else:
            self.export_fragment = export_fragment
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

        #  List of Database Fragments:   
        if self.list_db:
            self.list_dbentrys()
        ## Export !all! fragments   
        if self.export_all:
            self.export_all_fragments()
        ## Export one fragment         
        if self.export_fragment:
            self.do_export_fragment()
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
        export = Export(self.export_fragment)
        export.write_file()
        sys.exit(1)
    
    
    def list_dbentrys(self):
        '''
        list all entries in the db.
        '''
        gdb = global_DB()
        db = gdb.build_db_dict()
        print('\n Entries found in the databases:\n')
        print(' Fragment         | Line | DB Name    | Comment ')
        print(' ---------------------------------------------------------------------------')
    
        frags = sorted(db.keys())
        for i in frags:
            line = ' {:<17}| {:<5}| {:<11}| {}'.format(
                    i, gdb.get_line_number_from_fragment(i), 
                    gdb.get_db_from_fragment(i), 
                    gdb.get_comment_from_fragment(i))
            print(line[:79])
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
        mog = ImportGRADE(self.import_grade)
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
        gdb = global_DB()
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
        #print target_lines, res_target_atoms
        for i in target_lines:
            i = int(i)
            #reslist[i]=' '+reslist[i]
            rle.remove_line(i, rem=False, remove=False, frontspace=True)
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

    
    def chunks(self, l, n):
        """ Yield successive n-sized chunks from l.
        """
        for i in range(0, len(l), n):
            fehlt = 4 - len(l[i:i+n])
            i=i-fehlt
            yield l[i:i+n]

    
    def generate_flat_restraints(self, am, coords, cell):
        import networkx as nx
        from misc import vol_tetrahedron
        G = am.get_adjmatrix
        l = nx.cycle_basis(G)
        
        if not l:
            return False
        for ring in l:
            x=len(ring)
 #           print( ((x+3)//4)+(x%4)//4 , 'anzahl')
            parts = self.chunks(ring, 4)
            for i in parts:
                for fl in i:
                    print(coords[fl])
#            for n, atom in enumerate(ring):
#                atcoord = coords[atom]
#                x = atcoord
#                print(x, atom)
#                if n == 3:
#                    print('\n')
#                    break
#            #print(ring)
#        rem rem dsr put naphthalene with c8 c6 c7 on q6 Q4 q7 resi PART 3 occ -31 =
#         rem dfix
        #vol_tetrahedron(a, b, c, d, cell)
        #print(l)
        
    
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
        dfixes = re.get_formated_12_dfixes+re.get_formated_13_dfixes
        print(self.generate_flat_restraints(am, coords, cell))
        return ''.join(dfixes)
    
    
    def main(self): 
        '''
        main object to run DSR as command line program
        '''
        print(progname) 
        # The database content:
        gdb = global_DB()
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
        
        print('Inserting {} into res File.'.format(fragment))
        db_source_atoms = dsr_dict['source']
        print('Source atoms:', ', '.join(db_source_atoms))
        res_target_atoms = dsr_dict['target']
        print('Target atoms:', ', '.join(res_target_atoms))
        
        # several checks if the atoms in the dsr command line are consistent
        check_source_target(db_source_atoms, res_target_atoms, dbatoms)
        
        num = NumberScheme(reslist, dbatoms, resi.get_resinumber)
        numberscheme = num.get_fragment_number_scheme()
        afix = InsertAfix(reslist, dbatoms, dbtypes, dbhead, dsr_dict, 
                          sfac_table, find_atoms, numberscheme)
        afix_entry = afix.build_afix_entry()
        # line where the dsr command is found in the resfile:
        dsrline = dsrp.find_dsr_command(reslist) 
        reslist[dsrline] = reslist[dsrline]+'\n'
        reslist.insert(dsrline+1, afix_entry)
    
        ##### comment out all target atom lines in replace mode:  
        if dsr_dict['command'] == 'REPLACE':
            replacemode(res_target_atoms, rle, reslist)

        # write to file:
        basefilename = filename_wo_ending(self.res_file)
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
        
        if dsr_dict['dfix']:
            resinumber = resi.get_resinumber
            dfix = self.generate_dfix_restraints(lf, 
                                            reslist, 
                                            numberscheme, 
                                            resinumber,
                                            cell,
                                            dsr_dict['part'])
            if resinumber:
                for n, line in enumerate(reslist):
                    if line.upper().startswith('RESI'):
                        if line.split()[1] == str(resinumber):
                            line = '{} \n{}'.format(line, dfix) # insert restraints after residue definition
                            reslist[n] = line
            else:
                pass
                line = '{}'.format(dfix) # insert restraints after dsrline
                reslist[dsrline] = reslist[dsrline]+line
                
        if not self.no_refine:
            rl.write_resfile(reslist, '.res')
        
  
    
if __name__ == '__main__':
    '''main function'''
    dsr = DSR(no_refine=True)
    #dsr = DSR()

    

