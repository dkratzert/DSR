#-*- encoding: utf-8 -*-
#mÃ¶p
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
from atomhandling import Elem_2_Sfac, rename_restraints_atoms
import misc
import os
import sys
import constants
import fnmatch
from constants import RESTRAINT_CARDS
from misc import id_generator


__metaclass__ = type  # use new-style classes


def add_residue_to_dfix(dfix_head, resinum):
    """
    Add a residue to a list of DFIX/DANG restraints
    DFIX 1.234 C1 C2 -> DFIX 1.234 C1_4 C2_4
    >>> add_residue_to_dfix(['DFIX 1.456 C1 C2', 'DFIX 1.212 C3 C4'], 4)
    ['DFIX  1.456  C1_4  C2_4\\n', 'DFIX  1.212  C3_4  C4_4\\n']
    >>> add_residue_to_dfix(['DFIX 1.456 C1 C2', 'DFIX 1.212 C3 C4'], '5')
    ['DFIX  1.456  C1_5  C2_5\\n', 'DFIX  1.212  C3_5  C4_5\\n']
    >>> add_residue_to_dfix(['FLAT C6 C1 C2 C3', 'FLAT  C7  C1  C2  C3'], '2')
    ['FLAT  C6_2  C1_2  C2_2  C3_2\\n', 'FLAT  C7_2  C1_2  C2_2  C3_2\\n']
    """
    newhead = []
    for line in dfix_head:
        line = line.split()
        first = line[0]
        for num, item in enumerate(line):
            try:
                int(item[0])
            except:
                line[num] = line[num] + '_' + str(resinum)
                continue
        line = first+'  '+'  '.join(line[1:]) + '\n'
        newhead.append(line)
    return newhead


def insert_dsr_warning():
    """
    information to insert into .res-file
    TODO: insert this text only once in the file and also the
          comments from GRADE
    """
    txt = "rem the following was inserted by DSR:\n" \
          "http://dx.doi.org/10.1107/S1600576715005580\n"
    return txt


class Afix(object):
    """
    methods for the AFIX entry
    - dbhead is modified by Resi() if residues are used!
      RESI num class ist inserted there
    """

    def __init__(self, reslist, dbatoms, fragment_atom_types, restraints, dsrp, sfac_table,
                 find_atoms, numberscheme, options, dfix_restraints=False):
        """
        :param reslist:      list of the .res file
        :type reslist: list
        :param dbatoms:      list of the atoms in the database entry
                             [['O1', 3, '-0.01453', '1.66590', '0.10966'], ... ]
        :type dbatoms: list
        :param fragment_atom_types:      ['N', 'C', 'C', 'C']
        :param restraints:       restraints from the database
        :param dsrp: dsr parser object
        :param sfac_table:   SHELXL SFAC table as list like: ['C', 'H', 'O', 'F', 'Al', 'Ga']
        :param find_atoms:   FindAtoms() object
        :param numberscheme: atoms numbering scheme like: ['O1A', 'C1A', 'C2A', 'F1A', 'F2A', 'F3A', 'C3A']

        """
        self._reslist = reslist
        self._find_atoms = find_atoms
        self._dbatoms = dbatoms
        self._restraints = restraints
        self.dfix_restraints = dfix_restraints
        self._fragment_atom_types = fragment_atom_types
        self._sfac_table = sfac_table
        self.numberscheme = numberscheme
        self.dsrp = dsrp
        self.part = dsrp.part
        self.occ = dsrp.occupancy
        self.source_atoms = dsrp.source
        self.target_atoms = dsrp.target
        self._dfix = dsrp.dfix
        self.options = options
        self.rand_id_dfx = id_generator(size=7)
        self.rand_id_afix = id_generator(size=7)

    def collect_all_restraints(self):
        """
        :return all_restraints: list
        collects all restraints in the resfile and returns a list with them
        [['RIGU_CF3', 'O1', '>', 'F9'], '...']
        """
        all_restraints =[]
        for n, resline in enumerate(self._reslist):
            resline = resline.strip(' \n\r')
            #resline = resline.split()
            try:
                resline[:4]
            except:
                continue
            if resline[:4] in RESTRAINT_CARDS:
                # see for the next  lines if the lines continues with "=":
                line = 0
                while resline[-1] == '=':
                    resline = resline[:-1]+self._reslist[n+line+1]
                    line += 1
                    if not resline[-1] == '=':
                        break
                    if line > 500:
                        break
                all_restraints.append(resline)
        return all_restraints

    @staticmethod
    def remove_duplicate_restraints(restraints, all_restraints, residue_class=''):
        """
        removes restraints from the header which are already
        in the res-file.

        :param restraints:         database header (list of strings)
        :param residue_class:  SHELXL residue class
        :type residue_class:   string
        :param all_restraints: all restraints in the res file
        :type all_restraints:  list
        :return new_restr: list
        >>> dbhead = ["SADI 0.02 C1 C2 C2 C3 C3 C4", "SADI 0.04 C1 C3 C3 C5", "DFIX 1.45 C1 C2"]
        >>> all_restraints = ["SADI C1 C2 C2 C3 C3 C4", "SADI C1 C3 C3 C5", "DFIX C1 C2", "SADI C4 C5 C5 C6"]
        >>> Afix.remove_duplicate_restraints(dbhead, all_restraints)
        <BLANKLINE>
        Already existing restraints were not applied again.
        ['', '', '']

        >>> all_restraints = ["SADI 0.02 C1 C2 C2 C3 C3 C4", "SADI 0.04 C1 C3 C3 C5", "DFIX 1.45 C1 C2"]
        >>> dbhead = ["SADI C1 C2 C2 C3 C3 C4", "SADI C1 C3 C3 C5", "DFIX C1 C2", "SADI C4 C5 C5 C6"]
        >>> Afix.remove_duplicate_restraints(dbhead, all_restraints)        
        <BLANKLINE>
        Already existing restraints were not applied again.
        ['', '', '', 'SADI C4 C5 C5 C6']
        """
        modified = False
        new_restr = restraints[:]
        for num, line in enumerate(restraints):
            line = Afix.remove_stddev_from_restraint(line.split())
            for restr in all_restraints:
                restr = restr.split()
                restr = Afix.remove_stddev_from_restraint(restr)
                if line == restr:
                    new_restr[num] = ''
                    modified = True
                    break
        if modified:
            if residue_class:
                print('\nAlready existing restraints for residue "{}" were not '
                      'applied again.'.format(residue_class))
            else:
                print('\nAlready existing restraints were not applied again.')
        return new_restr

    @staticmethod
    def remove_stddev_from_restraint(restr):
        # type: (list) -> list
        """
        Parameters
        ----------
        restr: list of restraints

        Returns list of restraints without standars deviation
        -------
        
        >>> r = ['SADI', '0.02', 'C1', 'C2', 'C3', 'C4']
        >>> r2 = ['SADI', 'C1', 'C2', 'C3', 'C4']
        >>> Afix.remove_stddev_from_restraint(r)
        ['SADI', 'C1', 'C2', 'C3', 'C4']
        >>> Afix.remove_stddev_from_restraint(r2)
        ['SADI', 'C1', 'C2', 'C3', 'C4']
        """
        new = []
        # find out where the atoms begin (leave out numbers):
        for num, i in enumerate(restr):
            if i[0].isalpha():
                new.append(i)
        return new

    @staticmethod
    def distance_and_other_restraints(restraints):
        """
        Devides restraints in distance restraints (distance)
        and all other lines (others)
        Restraints are instead inserted after fragment fit

        :param restraints:  database header
        """
        distance = []
        others = []
        for num, headline in enumerate(restraints):  # @UnusedVariable
            headline = headline.strip().split()
            try:
                headline[0]
            except IndexError:
                continue
            if headline[0][:4] in constants.DIST_RESTRAINT_CARDS:
                distance.append(' '.join(headline)+'\n')
            else:
                others.append(' '.join(headline)+'\n')
        return [distance, others]

    def combine_names_and_coordinates(self):
        """
        Combines the target atom names with the coordinates from the -target option.
        Douplicate q-peak names are explicitely allowed here for special positions.
        Therefore, the target atoms are named DUM0, DUM1, DUM2, ...
        :rtype: dict
        """
        atoms = {}
        chunk = misc.chunks(self.options.target_coords, 3)
        if len(chunk) != len(self.target_atoms):
            print("*** Different number of target atoms and target coordinates! Can not proceed. ***")
            sys.exit()
        tmp = self.target_atoms
        self.target_atoms = []
        for num, at in enumerate(tmp):
            # Handles douplicate q-peak names:
            at = "DUM{}".format(num)
            atoms[at] = chunk[num]
            self.target_atoms.append(at)
        return atoms

    def build_afix_entry(self, external_restraints, dfx_file_name, resi):
        """
        build an afix entry with atom coordinates from the target atoms

        :type resi: Resi
        :param external_restraints:  True/False decision if restraints should be
                                     written to external file
        :param dfx_file_name:        name of file for external restraints
        :param resi:                Residue() object
        """
        old_atoms = []
        afix_list = []   # the final list with atoms, sfac and coordinates
        e2s = Elem_2_Sfac(self._sfac_table)
        new_atomnames = list(reversed(self.numberscheme)) # i reverse it to pop() later
        # Residue is active:
        if self.dsrp.resiflag:
            self._restraints = resi.format_restraints(self._restraints)
            self._restraints = self.remove_duplicate_restraints(self._restraints, self.collect_all_restraints(),
                                                                resi.get_residue_class)
            if not external_restraints:
                self._restraints += ['RESI {} {}'.format(
                    resi.get_resinumber, resi.get_residue_class)]
        # No residue:
        else:
            # applies new naming scheme to head:
            old_atoms = [ i[0] for i in self._dbatoms]
            self._restraints = rename_restraints_atoms(new_atomnames, old_atoms, self._restraints)
            self._restraints = self.remove_duplicate_restraints(self._restraints, self.collect_all_restraints())
        # decide if restraints to external file or internal:
        distance_and_other = self.distance_and_other_restraints(self._restraints)
        distance = distance_and_other[0]
        other_head = distance_and_other[1]
        # External restraints:
        if external_restraints and not self.options.rigid_group:
            # in case of dfix, write restraints to file after fragment fit
            self._restraints = misc.wrap_headlines(distance)
            # returns the real name of the restraints file:
            if self.dfix_restraints:
                # DFIX enabled:
                pname = os.path.splitext(dfx_file_name)
                dfx_file_name = pname[0]+"_dfx"+pname[1]
                if self.dsrp.resiflag:
                    # External, no dfix but with residue:
                    self.dfix_restraints = add_residue_to_dfix(self.dfix_restraints, resi.get_resinumber)
                else:
                    # No residue but dfix and external:
                    self.dfix_restraints = rename_restraints_atoms(new_atomnames, old_atoms, self.dfix_restraints)
                dfx_file_name = self.write_dbhead_to_file(dfx_file_name, self.dfix_restraints, resi.get_residue_class,
                                                          resi.get_resinumber)
            else:
                # DFIX disabled:
                dfx_file_name = self.write_dbhead_to_file(dfx_file_name, self._restraints, resi.get_residue_class,
                                                          resi.get_resinumber)
                self._restraints = self._restraints = other_head
            if self.dfix_restraints:
                self._restraints = other_head
            if self.dsrp.resiflag:
                self._restraints += ['RESI {} {}'.format(resi.get_resinumber, resi.get_residue_class)]
        # No external restraints:
        else:
            if self.dfix_restraints:
                if self.dsrp.resiflag:
                    self.dfix_restraints = add_residue_to_dfix(self.dfix_restraints, resi.get_resinumber)
                self._restraints = other_head + self.dfix_restraints
                if not self.dsrp.resiflag:
                    self._restraints = rename_restraints_atoms(new_atomnames, old_atoms, self._restraints)
            self._restraints = misc.wrap_headlines(self._restraints)
        # list of atom types in reverse order
        reversed_fragm_atom_types = list(reversed(self._fragment_atom_types))
        if self.options.target_coords:
            # {'C1': ['1.123', '0.7456', '3.245']}
            coordinates = self.combine_names_and_coordinates()
        else:
            coordinates = self._find_atoms.get_atomcoordinates(self.target_atoms)
        occ_part = self.occ
        if not self.occ:
            occ_part = ''
            occ = "11.00"
        else:
            occ = self.occ
        # a list of zeroed atom coordinates (afix_list) is built:
        for i in self._dbatoms:
            line = []
            sfac_num = str(e2s.elem_2_sfac(reversed_fragm_atom_types.pop()))
            line.insert(0, str(i[0]))  # Atom name
            line.insert(1, sfac_num)   # SFAC number
            line.insert(2, '0       ')
            line.insert(3, '0       ')
            line.insert(4, '0       ')
            line.insert(5, occ)
            line.insert(6, '0.04')
            afix_list.append(line)
        # for every atoms both in afix list and source_atoms, change the
        # coordinates to the respective value of the target atom:
        ind = 0
        for n, i in enumerate(afix_list):
            if i[0].upper() in self.source_atoms:
                ind = self.source_atoms.index(i[0].upper())
                try:
                    afix_list[n][2:5] = coordinates[self.target_atoms[ind]]
                except IndexError:
                    print('*** More source than target atoms present! Exiting... ***')
                    sys.exit(False)
        new_atom_list = []
        new_atomnames.reverse()
        for n, i in enumerate(afix_list):
            i[0] = new_atomnames[n]
            i[2] = '{:>10.6f}'.format(float(i[2]))
            i[3] = '{:>10.6f}'.format(float(i[3]))
            i[4] = '{:>10.6f}'.format(float(i[4]))
            new_atom_list.append('   '.join(i).rstrip())
        atoms = '\n'.join(new_atom_list)
        afixnumber = '179'   # makes afix 179 default
        if self.part:
            part = 'PART '+str(self.part)+' '+str(occ_part)
            part2 = 'PART 0'
        else:
            part = ''
            part2 = ''
        if self.dsrp.resiflag:
            resi_end = 'RESI 0'
        else:
            resi_end = ''
        if external_restraints and not self.options.rigid_group:
            if self.dsrp.resiflag:
                self._restraints += '\nREM The restraints for residue {} are in this file:\nrem +{}\nREM {}\n'\
                                            .format(resi.get_residue_class, dfx_file_name, self.rand_id_dfx)
            else:
                self._restraints += '\nREM The restraints for this moiety are in this file:\nrem +{}\nREM {}\n'\
                                            .format(dfx_file_name, self.rand_id_dfx)
        if self.options.rigid_group:
            afixtag = ''
            if self.dsrp.resiflag:
                self._restraints = 'RESI {} {}\n'.format(resi.get_residue_class, resi.get_resinumber)
            else:
                self._restraints = ''
        else:
            afixtag = 'REM '+self.rand_id_afix
            self._restraints = ''.join(self._restraints)
        # Adds a "SAME_resiclass firstatom > lastatom" to the afix:
        same = ''
        if self.dsrp.resiflag and not self.dsrp.dfix and not self.options.rigid_group:
            same = "SAME_{} {} {}\n".format(resi.get_residue_class, new_atomnames[0], new_atomnames[-1])
        afix = '{0}{8}{1}\n' \
               'AFIX {2}\n' \
               '{3}\n' \
               '{4}\n' \
               'AFIX 0\n' \
               '{5}\n' \
               '{6}\n' \
               '{7}\n\n'.format(self._restraints,  # 0
                                part,  # 1
                                str(afixnumber),  # 2
                                afixtag,  # 3
                                atoms,  # 4
                                afixtag,  # 5
                                part2,  # 6
                                resi_end,  # 7
                                same)        # 8
        return afix

    def write_dbhead_to_file(self, filename, dbhead, resi_class, resi_number):
        """
        write the restraints to an external file
        :param filename:     filename of database file
        :param dbhead:       database header
        :param resi_class:   SHELXL residue class
        :param resi_number:  SHELXL residue number
        :return filename:    full file name where restraints will be written
        """
        number = '1'
        files = []
        # find a unique number for the restraint file:
        for filen in misc.sortedlistdir('.'):
            if fnmatch.fnmatch(filen, 'dsr_*_' + filename):
                filenum = filen.split('_')
                if str.isdigit(filenum[1]):
                    files.append(filenum[1])
        filepath, filename = os.path.split(os.path.abspath(filename))
        try:
            number = str(int(files[-1]) + 1)
        except IndexError:
            pass
        if not self.dsrp.resiflag:  # no residues
            filename = 'dsr_' + number + '_' + filename
        if self.dsrp.resiflag and resi_number:  # only residue number known
            filename = 'dsr_' + resi_class + '_' + resi_number + '_' + filename
        if self.dsrp.resiflag and not resi_number and resi_class:  # only residue class known
            filename = 'dsr_' + resi_class + '_' + filename
        if os.path.isfile(os.path.abspath(filename)):
            print('Previous restraint file found. Using restraints from "{}"'.format(filename))
            return filename
        try:
            dfix_file = open(os.path.join(filepath, filename), 'w')  # open the ins file
        except IOError:
            print('*** Unable to write restraints file! Check directory write permissions. ***')
            sys.exit(False)
        print('Restraints were written to "{}"'.format(os.path.join(filepath, filename)))
        for i in dbhead:  # modified reslist
            dfix_file.write("%s" % i)  # write the new file
        dfix_file.close()
        return filename

if __name__ == '__main__':
    import doctest
    failed, attempted = doctest.testmod()#verbose=True)
    if failed == 0:
        print('passed all {} tests!'.format(attempted))







