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
import pstats
from dsrparse import DSR_Parser
from options import OptionsParser
from dbfile import global_DB, ImportGRADE
from resfile import ResList, ResListEdit
from atomhandling import SfacTable, get_atomtypes, check_source_target, Elem_2_Sfac
from atomhandling import FindAtoms
from afix import InsertAfix
from refine import ShelxlRefine
from resi import Resi
from misc import get_replace_mode, find_line
from restraints import ListFile, Lst_Deviations


# TODO and ideas:
# -To the manual: Avogradro/res-file -> rename -> mercury -> mol2-file -> GRADE
# -To manual: what to do if something does not work.
# -more detailed comments for grade imports:
#   REM Produced by Grade Web Server http://grade.globalphasing.org
#   REM GEN: Generated by GRADE 1.2.5 (December 20 2013)
#   REM GEN: from SMILES C1=CC=C2C(=C1)C=C(C3=CC=CC=C23)Br
#   REM GEN: using MOGUL 1.6(RC5), CSD as535be with quantum mechanics RM1
#   REM grade-cif2shelx output
#   REM Version: 0.0.5 <Dec 20 2013>
# -make a proper python module from DSR
# -try to iterate all target atoms if fit is really bad
# -Iteratively analyze disagreeable restraints and set ESDs of relative
#  restraints accordingly. (Only if shift is not too big, has to be nearly converged)
# -select atoms in shelxle > export to db
# -if part > 2 and same residue and positive occ then SUMP
# -use +filename.dfix for restraints.
# -If second part of a disordered atom is very close to the first, make EADP?
# -If e.g. OCC in dsrline and index(OCC)+1 != ' ', then dsrline.insert(' ', index(OCC)+1)
# -check for residues with same class and differing atom names inside.
# -should I make rem dsr IMPORT from Atom1 to Atom2 ?
#  with or without their restraints? RESI class might work?
# -add SIMU and RIGU after Grade import
# -debian package: /usr/src/packages/BUILD # dpkg-deb --build dsr


VERSION = '1.2.12'
progname = '\n----------------------------- D S R - v{} ----------------------------------'.format(VERSION)

def export_fragment(options):
    ''' 
    Exports the current fragment.
    '''
    from export import Export
    export = Export(options)
    export.write_file()
    sys.exit(1)


def list_dbentrys():
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
    print('\n Feel free to add more fragments to "dsr_user_db.txt" in the program directory\n or mail them to dkratzert@gmx.de.\n')
    
    for fragment in list(db.keys()):
        gdb.check_consistency(db[fragment], fragment)
        gdb.check_db_atom_consistency(db[fragment]['atoms'], fragment)
        gdb.check_db_header_consistency(db[fragment]['head'], fragment)
    sys.exit()


def import_from_grade(options):
    '''
    imports a fragment from the GRADE webserver.
    '''
    mog = ImportGRADE(options.import_grade)
    mog.write_user_database()
    sys.exit(1)


def set_post_refine_cycles(shx, cycles):
    '''
    set the number of refinement cycles
    cycles must be a string
    '''
    try:
        shx.set_refinement_cycles(cycles)
    except(IndexError):
        print('Unable to set refinement cycles')
    shx.remove_afix()   # removes the afix 9
    

def export_all_fragments(options):
    '''
    export all database entries at once
    '''
    from export import Export
    gdb = global_DB()
    db = gdb.build_db_dict()
    dbnames = list(db.keys())
    for name in dbnames:
        export = Export(options, fragment=name)
        export.write_file()
    sys.exit(1)
    


def set_final_db_sfac_types(dbtypes, dbatoms, sfac_table):
    '''
    corrects the sfac types of the dbentry according to sfac card of the
    res file
    '''
    e2s = Elem_2_Sfac(sfac_table)
    atype = list(reversed(dbtypes))
    for line in dbatoms:                             # go through db entry
        line[1] = e2s.elem_2_sfac(atype.pop())       # replace scattering factor (line[1]) with true one


def replacemode():
    '''
    Target atoms are being replaced if this is executed
    '''
    print('Replace mode active.')
    target_lines = fa.get_atom_line_numbers(res_target_atoms)
    #print target_lines, res_target_atoms
    for i in target_lines:
        i = int(i)
        rle.remove_line(i, rem=False, remove=False, frontspace=True)
    fa.remove_adjacent_hydrogens(res_target_atoms)    


def go_refine(shx):
    '''
    actually starts the fragment fit
    '''
    try:
        shx.run_shelxl()
    except() as e:
        print(e)
        sys.exit()


def main(): 
    # options from the commandline options parser
    options = OptionsParser(progname)
    # The database content:
    gdb = global_DB()
    
    print(progname) # prints the version string on screen
    
    #  List of Database Fragments:   
    if options.list_db:
        list_dbentrys()
   
    ## Export !all! fragments   
    if options.export_all:
        export_all_fragments(options)

    ## Export one fragment         
    if options.export_fragment:
        export_fragment(options)
 
    ## Import a GRADE fragment          
    if options.import_grade:
        import_from_grade(options)

    rl = ResList(options.res_file)
    reslist = rl.get_res_list()
    find_atoms = FindAtoms(reslist)
    rle = ResListEdit(reslist, find_atoms)
    dsrp = DSR_Parser(reslist, rle)
    dsr_dict = dsrp.parse_dsr_line()
    fvarlines = rle.find_fvarlines()
    
    if dsrp.occupancy:
        rle.set_fvar(dsrp.occupancy, fvarlines)
    
    fragment = dsrp.fragment

    # line of the fragment card as string:
    fragline = gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
    dbatoms = gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
    dbhead = gdb.get_head_from_fragment(fragment)        # this is only executed once
    residue = gdb.get_resi_from_fragment(fragment)
    dbtypes = get_atomtypes(dbatoms)                 # the atomtypes of the dbentry as list e.g. ['C', 'N', ...]
    resi = Resi(reslist, dsr_dict, dbhead, residue, find_atoms)
    dbhead = resi.make_resihead()
    sf = SfacTable(reslist, dbtypes)
    sfac_table = sf.set_sfac_table()                 # from now on this sfac table is set
    
    ### corrects the atom type according to the previous defined global sfac table:
    set_final_db_sfac_types(dbtypes, dbatoms, sfac_table)


    ## Insert FRAG ... FEND entry: 
    rle.insert_frag_fend_entry(dbatoms, fragline, fvarlines)
    
  #  ######################################################
  #  ##    Output debug infos                            ##
  #  ######################################################   
  #
  #  if options.debug:
  #      try:
  #          p = pstats.Stats('dsr.profile')
  #      except(IOError):
  #          print '\nunable to find "dsr.profile"'
  #          sys.exit(-1)
  #      p.strip_dirs().sort_stats('cumulative').print_stats(12)
  #      p.strip_dirs().sort_stats('time').print_stats(12)
  #      sys.exit(1)
  #  
  #  #######################################################

    print('Inserting {} into res File.'.format(fragment))
    db_source_atoms = dsr_dict.get('source')
    print('Source atoms:', ', '.join(db_source_atoms))
    res_target_atoms = dsr_dict.get('target')
    print('Target atoms:', ', '.join(res_target_atoms))
    
    
    # several checks if the atoms in the dsr command line are consistent
    check_source_target(db_source_atoms, res_target_atoms, dbatoms)
    

    afix = InsertAfix(reslist, dbatoms, dbtypes, dbhead, dsr_dict, sfac_table, find_atoms)
    afix_entry = afix.build_afix_entry()
    # line where the dsr command is found in the resfile:
    dsrline = dsrp.find_dsr_command(reslist) 
    #reslist[dsrline-1] = reslist[dsrline-1]+'\n'
    reslist[dsrline] = reslist[dsrline]+'\n'
    reslist.insert(dsrline+1, afix_entry)

    ##### comment out all target atom lines in replace mode:  
    if dsr_dict.get('command') == 'REPLACE':
        replacemode()
    

    #############################################################
    ##    write to file:                                       ##
    basefilename = rl.filename_wo_ending(options.res_file)
    shx = ShelxlRefine(reslist, basefilename, find_atoms)
    if not options.no_refine:
        shx.set_refinement_cycles('0')
    shx.remove_acta()
    
    rl.write_resfile(reslist, '.ins')  
    
    
    ###########################################################
    ###  Refine with L.S. 0 to insert the fragment          ###
    if not options.no_refine:
        go_refine(shx)
    
    ### Display the results from the list file:
    lf = ListFile()
    lst_file = lf.read_lst_file()
    lfd = Lst_Deviations(lst_file)
    lfd.print_deviations()
    
    ### open res file again to restore 10 refinement cycles
    rl = ResList(options.res_file)
    reslist = rl.get_res_list()
    shx = ShelxlRefine(reslist, basefilename, find_atoms)
    if not options.no_refine:
        set_post_refine_cycles(shx, '8')
    
    if not options.no_refine:
        rl.write_resfile(reslist, '.res')
    





if __name__ == '__main__':
    '''main function'''
    
    main()
    

