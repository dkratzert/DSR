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
# each residue number must always have the same class!
# e.g 1 can be only TOL, not also NET4 but Tol can be 1 or 2.
#
# Wenn das dict leer ist werden überhaupt keine RESI betreffenden dinge in das insfile übernommen

from __future__ import print_function

import sys

from src.atomhandling import FindAtoms
from src.constants import RESTRAINT_CARDS
from src.dsrparse import DSRParser


def remove_resi(head):
    # type: (list) -> list
    """
    removes all resi commands and classes from restraints
    :param head: database header of the current fragment
    :type head: list
    """
    rhead = []  # head without resi
    delhead = []
    for dummy, line in enumerate(head):
        line = line.split()
        try:
            if line[0][:4] in RESTRAINT_CARDS:
                line[0] = line[0].split('_')[0]
        except IndexError:
            continue
        line = ' '.join(line)
        delhead.append(line)
    for line in delhead:
        line = line.strip(' \n\r').upper()
        if line.startswith('RESI'):
            continue
        rhead.append(line)
    return rhead


class Resi(object):
    def __init__(self, dsrp, db_residue_string, find_atoms):
        # type: (DSRParser, str, FindAtoms) -> NotImplemented
        """
        Handles the RESI instructions and restraints
        :param dsrp: Dsrparse object
        :type dsrp: object
        :param reslist: res file as list
        :type reslist: list
        :param dbhead: database header with restraints and residue
        :type dbhead: list
        :param db_residue_string: db residue string
        :type db_residue_string: string
        :param find_atoms: find_atoms object
        :type find_atoms: object
        """
        self._atoms_in_reslist = find_atoms.atoms_as_residues()
        self._residues_in_res = sorted(list(self._atoms_in_reslist))
        self.dsrp = dsrp
        self._resi_dict_dsr_command = {'class': None, 'number': None, 'alias': None}
        if self.dsrp.resiflag:
            # use a residue
            if self.dsrp.resi:
                # a residue was defined (class, number, or both)
                self._resi_dict_dsr_command = self.get_resi_syntax(self.dsrp.resi)
            else:
                # no residue defined, only active (use class from DB)
                self._resi_dict_dsr_command = {'class': None, 'number': None, 'alias': None}
        try:
            # use residue class from db if not given in dsr command:
            self._db_resi_list = db_residue_string.split()
        except AttributeError:
            print('No valid residue "RESI classname" in the database entry '
                  'of {} found.'.format(self.dsrp.fragment))
            sys.exit()
        self._resi_dict_db = self.get_resi_syntax(self._db_resi_list)
        self._combined_resi = self.build_up_residue()

    @property
    def get_resinumber(self):
        """
        Returns the residue number of the currently fitted fragment
        :rtype self._combined_resi['number']: string
        """
        return self._combined_resi['number']

    @property
    def get_residue_class(self):
        """
        Returns the residue class of the currently fitted fragment. Also is an indicator 
        if residues are active.
        :rtype self._combined_resi['class']: str
        """
        return self._combined_resi['class']

    def format_restraints(self, head):
        """
        in case of RESI, format the restraints like "SAME_class"
        """
        newhead = []
        for line in head:
            line = line.upper()
            line = line.split()
            try:
                line[0]
            except:
                continue
            if line[0] in RESTRAINT_CARDS:
                line[0] = line[0] + '_' + self._combined_resi['class']
                line = ' '.join(line)
            else:
                line = ' '.join(line)
            newhead.append(line)
        return newhead

    def get_unique_resinumber(self, resinum):
        # type: (str) -> str
        """
        Finds a unique resi number. If the number is already unique
        the given is used.
        :param resinum: residue number of the fragment
        :type resinum: string
        :return resinum: unique residue number
        """
        new_num = '1'
        if resinum == '0':
            while new_num in self._residues_in_res:
                new_num = str(int(new_num) + 1)
            return new_num
        elif resinum in self._residues_in_res and not resinum == '0':
            print('Warning: The residue number "{}" you have chosen is already ' \
                  'in use!'.format(resinum))
            while new_num in self._residues_in_res:
                new_num = str(int(new_num) + 1)
            return new_num
        else:
            return resinum

    def build_up_residue(self):
        # type: () -> dict
        """
        Decides which class and residue number should be used for the fragment.
        Returns a final dict with the residue settings.
        self._resi_dict_dsr_command is False if resi is enabled but no values given.
        """
        final_residue = {'class': None, 'number': None, 'alias': None}
        resiclass = None
        resinum = '0'
        resialias = None

        #### for the db entry
        if self._resi_dict_db['alias']:
            resialias = self._resi_dict_db['alias']
        if self._resi_dict_db['class']:
            resiclass = self._resi_dict_db['class'].upper()
        if self._resi_dict_db['number']:
            resinum = self._resi_dict_db['number']
        #### for the comlist entry
        if not self._resi_dict_dsr_command:
            return final_residue
        else:
            if self._resi_dict_dsr_command['alias']:
                resialias = self._resi_dict_dsr_command['alias']
            if self._resi_dict_dsr_command['class']:
                resiclass = self._resi_dict_dsr_command['class'].upper()
            if self._resi_dict_dsr_command['number']:
                resinum = self._resi_dict_dsr_command['number']
        final_residue = {'class': resiclass, 'number': resinum, 'alias': resialias}
        if resinum == '0':
            final_residue['number'] = self.get_unique_resinumber(resinum)
            print('No residue number was given. Using residue number {}.'.format(final_residue['number']))
        else:
            final_residue['number'] = self.get_unique_resinumber(resinum)
            print('Using residue number {}.'.format(final_residue['number']))
        if final_residue['class']:
            return final_residue

    @staticmethod
    def _wrong_syntax():
        print('This is not a valid RESIdue syntax!')

    def get_resi_syntax(self, inresi):
        """
        Checks if resi class, number and/or alias are present and valid.
        start with digit-> rest auch digit-> resinumber or alias
        start with letter-> rest letter or digit -> residue class

        The return value of this method is {'class', 'number', 'alias'}
        as a dictionary. Alias is empty if not given.

        The return value of just "RESI" in the command line is an empty dict
        :param resi: residue definition like ['3', 'CF3']
        :type resi: list
        """
        resi = inresi[:]
        resi_dict = {
            'class': '',
            'number': '',
            'alias': ''}
        if not resi:
            print('No valid RESI instruction found in the database entry!')
            sys.exit()
        try:
            resi.sort()
        except AttributeError:
            return resi_dict
        # any character of class must be a letter (New syntax in SHELXL 2017/1):
        if any([x.isalpha() for x in resi[-1]]):
            resi_dict['class'] = resi.pop()
        if len(resi) > 0:
            if str.isdigit(resi[0][0]):  # first character of number must be an digit
                resi_dict['number'] = resi[0]
                del resi[0]
            else:
                self._wrong_syntax()
                # not totally correct, but less difficult:
                print('Only digits are allowed for the residue number.')
                sys.exit()
        if len(resi) > 0:
            if str.isdigit(resi[0]):
                resi_dict['alias'] = resi[0]
                del resi[0]
            else:
                self._wrong_syntax()
                print('Only digits are allowed in the residue alias.')
                sys.exit()
        if resi_dict['number']:
            if len(resi_dict['number']) > 4:
                self._wrong_syntax()
                print('Only four digits allowed in residue number!')
                sys.exit()
        return resi_dict


if __name__ == '__main__':
    pass
