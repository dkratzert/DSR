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
from atomhandling import *
import misc
import constants


class InsertAfix(object):
    '''
    methods for the AFIX entry
    - dbhead is modified by Resi() if residues are used! RESI num class ist inserted there
    '''
    
    def __init__(self, reslist, dbatoms, dbtypes, dbhead, dsr_line, sfac_table, find_atoms):
        self.__reslist = reslist
        self._find_atoms = find_atoms
        self.__dbatoms = dbatoms
        self.__dbhead = dbhead
        self.__dbtypes = dbtypes
        self.__sfac = sfac_table
        self.part = dsr_line['part']
        self.occ = dsr_line['occupancy']
        #self.afixnumber = dsr_line['afix']
        self.source_atoms = dsr_line['source']
        self.target_atoms = dsr_line['target']
        self.__resi     = dsr_line['resi']
    

    def insert_dsr_warning(self):
        warn = 'rem the following was inserted by DSR:\n'
        return warn
    

    def remove_duplicate_restraints(self, dbhead):
        '''
        removes restraints from the header which are already 
        in the res-file
        '''
        newhead = dbhead[:]
        for resline in self.__reslist:
            resline = resline.strip().split()
            for num, headline in enumerate(dbhead):
                headline = headline.strip().split()
                if headline == resline and headline[0][:4] in constants.RESTRAINT_CARDS:
                    newhead[num] = ''
                    break
        return newhead
        
    
    def build_afix_entry(self): 
        '''
        build an afix entry with coordinates from the targetatoms
        '''
        atype = []       # list of atomtypes in reverse order
        afix_list = []   # the final list with atoms, sfac and coordinates
        e2s = Elem_2_Sfac(self.__sfac)
        # Fixme: move this to _init_
        num = NumberScheme(self.__reslist, self.__dbatoms, self.__resi)
        real_atomnames = list(reversed(num.get_fragment_number_scheme())) # i reverse it to pop() later
        # all non-atoms between start tag and FRAG card with new names:
        dbhead = rename_dbhead_atoms(real_atomnames, self.__dbatoms, self.__dbhead)
        dbhead = self.remove_duplicate_restraints(dbhead)
        atype = list(reversed(self.__dbtypes))
        coord = self._find_atoms.get_atomcoordinates(self.target_atoms)
        target = self.target_atoms[:]  # a copy because we edit it later
        # a list of zeroed atomcoordianes (afix_list) is built:
        for i in self.__dbatoms:
            l = []
            l.insert(0, str(i[0]))         # Atomname
            l.insert(1, str(e2s.elem_2_sfac(atype.pop())))  # SFAC number
            l.insert(2, '0       ')
            l.insert(3, '0       ')
            l.insert(4, '0       ')
            l.insert(5, '11.000')
            l.insert(6, '0.04')
            afix_list.append(l)
        # for every atoms both in afix list and source_atoms, change the
        # coordinates to the respective value of the target atom:
        ind = 0
        for n, i in enumerate(afix_list):
            if i[0].upper() in self.source_atoms:
                ind = self.source_atoms.index(i[0].upper())
                afix_list[n][2:5] =  coord[self.target_atoms[ind]]
        newlist = []
        for i in afix_list:
            i[0] = real_atomnames.pop()
            i[2] = i[2].ljust(8, '0').rjust(9, ' ') 
            i[3] = i[3].ljust(8, '0').rjust(9, ' ')
            i[4] = i[4].ljust(8, '0').rjust(9, ' ')
            newlist.append('    '.join(i).rstrip())
            
        atoms = '\n'.join(newlist)

        #atoms = misc.ll_to_string(afix_list) # make list of lists to string
            
        #if not self.afixnumber:  
        self.afixnumber = '179'   # makes afix 179 default 
        
        if not self.occ:
            self.occ = ''
        if self.part:
            part = 'PART '+str(self.part)+' '+str(self.occ)+'\n'
            part2 = 'PART 0\n'
        else:
            part = ''
            part2 = ''
        if ''.join(self.__resi):
            resi = 'RESI 0\n'
        else:
            resi = ''
        dbhead = ''.join(dbhead)
        warn = self.insert_dsr_warning()
        afix = warn+dbhead+str(part)+'AFIX '+str(self.afixnumber)+(
                '\n'+atoms+'\nAFIX 0\n'+part2+resi+'rem End of DSR entry\n\n')
        return afix


if __name__ == '__main__':
    from dbfile import global_DB
    from dsrparse import DSR_Parser
    from atomhandling import SfacTable
    from resfile import ResList, ResListEdit
    from resi import Resi
    
    options = OptionsParser()
    
    rl = ResList(options.res_file)
    reslist = rl.get_res_list()
    dsrp = DSR_Parser(reslist, rl)
    dsr_dict = dsrp.parse_dsr_line()
    find_atoms = FindAtoms(reslist)
    rle = ResListEdit(reslist, find_atoms)
    gdb = global_DB()
    db = gdb.build_db_dict()
    fragment = 'oc(cf3)3'
    fragline = gdb.get_fragline_from_fragment(fragment)  # full string of FRAG line
    dbatoms = gdb.get_atoms_from_fragment(fragment)      # only the atoms of the dbentry as list
    dbhead = gdb.get_head_from_fragment(fragment)        # this is only executed once
    residue = gdb.get_resi_from_fragment(fragment)
    dbtypes = get_atomtypes(dbatoms)
    resi = Resi(reslist, dsr_dict, dbhead, residue, find_atoms)
    dbhead = resi.make_resihead()

    
   
    sf = SfacTable(reslist, dbtypes)
    sfac_table = sf.set_sfac_table()
    
    afix = InsertAfix(reslist, dbatoms, dbtypes, dbhead, dsr_dict, sfac_table, find_atoms)
    print afix.build_afix_entry()




