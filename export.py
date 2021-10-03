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

from __future__ import print_function

from copy import deepcopy

from atomhandling import get_atomtypes
from atoms import Element
from dbfile import ParseDB, invert_atomic_coordinates
from misc import wrap_headlines, wrap_stringlist
from restraints import Restraints

__metaclass__ = type  # use new-style classes


class Export():
    """
    This class implements the export of a database entry to a .res file.
    Included are the minimal informations which are needed to get a valid res file.
    e.g.:
    TITL toluene
    CELL  0.71  11.430  12.082  15.500  106.613  100.313  90.68
    SFAC C
    UNIT 1
    WGHT 0.05
    FVAR 1.0
    C1   1  0.268330  0.478380  0.161680  11.00   0.04
    C2   1  0.205960  0.555770  0.217990  11.00   0.04
    C3   1  0.249400  0.600760  0.310040  11.00   0.04
    C4   1  0.357300  0.568990  0.348900  11.00   0.04
    C5   1  0.420800  0.492470  0.294060  11.00   0.04
    C6   1  0.376630  0.447580  0.201340  11.00   0.04
    C7   1  0.221500  0.430400  0.060360  11.00   0.04
    HKLF 4
    END
    """

    def __init__(self, gdb, invert=False):
        """
        :param fragment_name: string, name of the database fragment
        :param invert:        bool, should the coordinates be inverted?
        """
        self.invert = invert
        self._gdb = gdb

    def expand_calced_cell(self, cell, atoms):
        # type: (list, list) -> list
        """
        TODO: make this crude hack more elegant!
        In calculated structure the cell is 1 1 1 90 90 90. Shelxle has problems
        with that when growing. So the cell is expanded to 50 50 50
        >>> gdb = ParseDB('../dsr_db.txt')
        >>> exp = Export(gdb)
        >>> atoms = gdb.get_atoms('benzene')
        >>> exp.expand_calced_cell([1, 1, 1, 90, 90, 90], atoms)
        [[50, 50, 50, 90, 90, 90], [['C1', 1, 0.0176, -0.006618, 0.005344], ['C2', 1, 0.015724, -0.007554, 0.004762], ['C3', 1, 0.015212, -0.006368, 0.003851], ['C4', 1, 0.016584, -0.004244, 0.00351], ['C5', 1, 0.018464, -0.003288, 0.00408], ['C6', 1, 0.018976, -0.004464, 0.004998]]]
        """
        atoms = deepcopy(atoms)
        summe = int(sum(cell[0:3]))  # this is to detect calculated structures
        if summe == 3:  # 1+1+1=3!
            for coord in range(2, 5):  # x, y, z of coordinates
                for line in atoms:  # for every atom line
                    line[coord] = round(line[coord] / 50, 6)
            # now the new 50,50,50 cell:
            for n in range(0, 3):
                cell[n] = 50
        # cell = [str(x) for x in cell]
        return [cell, atoms]

    def make_dfix(self, fragname):
        fragname = fragname.lower()
        restr = Restraints(fragname, self._gdb)
        dfix_12 = restr.get_formated_12_dfixes()
        dfix_13 = restr.get_formated_13_dfixes()
        flats = restr.get_formated_flats()
        dfix_head = dfix_12 + dfix_13 + flats
        return dfix_head

    def export_all_fragments(self):
        """
        Export all database entries at once
        """
        import sys
        for fragment in self._gdb:
            self.write_res_file(fragment)
        sys.exit(1)

    def export_resfile(self, fragname):
        """
        exports a .res file from a database entry to be viewed in a GUI
        >>> db = ParseDB('../dsr_db.txt')
        >>> ex = Export(db)
        >>> ex.export_resfile('water') # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        ['TITL water\\n', 'REM This file was exported by DSR version ...\\n',
        'REM Name: Water, H2O\\n', 'REM Source: pbe1pbe/6-311++G(3df,3pd), Ilia A. Guzei\\n',
        '\\n', 'CELL 0.71073  50.0000  50.0000  50.0000  90.0000  90.0000  90.0000\\n',
        'ZERR    1.00   0.000    0.000    0.000    0.000    0.000    0.000\\n',
        'LATT  -1\\n', 'SFAC O  H\\n', 'UNIT 1  1 \\n', 'REM  RESIDUE: H2O\\n',
        'REM Sum formula: H2 O1 \\n', 'WGHT  0.1\\n', 'FVAR  1.0\\n',
        'rem Restraints from DSR database:\\n',
        'DFIX 0.9584 0.001 O1 H1 O1 H2\\nDFIX 1.5150 0.001 H1 H2\\n',
        'rem Restraints from atom connectivities:\\n',
        ['DFIX 0.9584 H2   O1  \\n', 'DFIX 0.9584 H1   O1  \\n',
        'DANG 1.5151 H1   H2  \\n'], 'rem end of restraints\\n', '\\n',
        ['O1   1     0.00000   0.00000   0.00000   11.0   0.04\\n',
        'H1   2     0.01917   0.00000   0.00000   11.0   0.04\\n',
        'H2   2    -0.00478   0.01856   0.00000   11.0   0.04\\n'],
        '\\nHKLF 0\\nEND\\n']
        """
        fragname = fragname.lower()
        cell = self._gdb.get_cell(fragname)
        atoms = self._gdb.get_atoms(fragname)
        # expands the cell of calculated structures:
        cell, atoms = self.expand_calced_cell(cell, atoms)
        cellstring = '{:>8.4f} {:>8.4f} {:>8.4f} {:>8.4f} {:>8.4f} {:>8.4f}'.format(*cell)
        if self.invert:
            print("Fragment inverted.")
            atoms = invert_atomic_coordinates(atoms)
        try:
            from dsr import VERSION
        except ImportError:
            VERSION = ''
        sfac = []
        res_export = []
        self._gdb.check_consistency(fragname)
        self._gdb.check_db_atom_consistency(fragname)
        for i in get_atomtypes(atoms):  # build sfac table from atomtypes
            if i not in sfac:
                sfac.append(i)
        atlist = []
        for i in get_atomtypes(atoms):  # atomtypes in the db_entry
            for y, x in enumerate(sfac):
                if x == i:
                    atlist.append(y + 1)
        for n, i in enumerate(atlist):
            atoms[n][1] = i
        # build the UNIT table:
        unit = []
        for i in sfac:
            unit.append('1 ')  # no matter what number
        final_atomlist = [('{:4.4s} {:4.2s} {:>8.5f}  {:>8.5f}  {:>8.5f}   11.0   0.04\n'.format(
            str(i[0]), str(i[1]), float(i[2]), float(i[3]), float(i[4]))) for i in atoms]
        res_export.append('TITL ' + fragname + '\n')  # title card with fragment name
        try:
            res_export.append('REM This file was exported by DSR version {}\n'.format(VERSION))
        except NameError:
            pass
        comments = self._gdb[fragname]['comments']
        comments = wrap_stringlist(comments, 75)
        source = self._gdb[fragname]['source']
        source = '\n'.join(wrap_stringlist([source], 70))
        name = self._gdb[fragname]['name']
        name = '\n'.join(wrap_stringlist([name], 75))
        res_export.append('REM Name: {}'.format(name))
        res_export.append('REM Source: {}'.format(source))
        res_export.append("{}\n".format('\n'.join(comments)))
        res_export.append('CELL 0.71073 ' + cellstring + '\n')  # the cell with wavelength
        res_export.append('ZERR    1.00   0.000    0.000    0.000    0.000    0.000    0.000\n')
        res_export.append('LATT  -1\n')
        res_export.append('SFAC ' + '  '.join(sfac) + '\n')
        res_export.append('UNIT ' + ' '.join(unit) + '\n')
        res_export.append('REM  RESIDUE: {}\n'.format(self._gdb.get_resi(fragname)))
        res_export.append('REM Sum formula: {}\n'.format(self._gdb.get_sum_formula(fragname)))
        res_export.append('WGHT  0.1' + '\n')
        res_export.append('FVAR  1.0' + '\n')
        try:
            res_export.append('rem Restraints from DSR database:\n')
            res_export.append(''.join(wrap_headlines(self._gdb.get_restraints(fragname))))
        except:
            pass
        try:
            res_export.append('rem Restraints from atom connectivities:\n')
            res_export.append(self.make_dfix(fragname))
            res_export.append('rem end of restraints\n')
        except Exception as e:
            print("*** Error during restraints generation: {} ***".format(e))
            pass
        res_export.append('\n')
        res_export.append(final_atomlist)  # the atoms
        res_export.append('\nHKLF 0\nEND\n')  # the end
        return res_export

    def copy_to_clipboard(self, fragname):
        # type: (str) -> bool
        """
        copys the exported atoms to the clipboard including FRAG  FEND commands

        Export example:

        FRAG
        C1   1     1.2000  -0.0230   3.6150
        C2   1     1.2030  -0.0120   2.1060
        C3   1     0.0150  -0.0110   1.3900
        C4   1     0.0150  -0.0010   0.0050
        C5   1     1.2080   0.0080  -0.6880
        C6   1     2.3980   0.0060   0.0090
        C7   1     2.3940  -0.0040   1.3940
        FEND
        """
        import pyperclip
        fragname = fragname.lower()
        clip_text = []
        cell = self._gdb.get_cell(fragname)
        atoms = self._gdb.get_atoms(fragname)
        atoms = self.format_atoms_for_export(cell, atoms, False)
        atoms = '\n'.join(atoms)
        clip_text.append('FRAG')
        clip_text.append('\n' + atoms)
        clip_text.append('\nFEND')
        text = ' '.join(clip_text)
        pyperclip.setcb(text)
        return True

    @staticmethod
    def format_atoms_for_export(cell, atoms, gui=False):
        # type: (list, list, bool) -> list
        """
        Returns properly formated cartesian coordinates for fragment export.
        Atom;;number;;x;;y;;z

        >>> gdb = ParseDB('../dsr_db.txt')
        >>> atoms = gdb.get_atoms('toluene')
        >>> cell = gdb.get_cell('toluene')
        >>> exp = Export(gdb=gdb)
        >>> exp.format_atoms_for_export(cell, atoms, gui=False) # doctest: +NORMALIZE_WHITESPACE +REPORT_NDIFF
        ['C1      6  1.7810   7.1491  12.0042',
         'C2      6  2.2009   8.3068  11.1376',
         'C3      6  1.2689   9.0217  10.3903',
         'C4      6  1.6422  10.0777   9.5884',
         'C5      6  2.9808  10.4443   9.5173',
         'C6      6  3.9205   9.7497  10.2541',
         'C7      6  3.5389   8.6909  11.0530']
        >>> exp.format_atoms_for_export(cell, atoms, gui=True) # doctest: +NORMALIZE_WHITESPACE +REPORT_NDIFF
        ['C1 6 1.78099 7.14907 12.00423', 'C2 6 2.20089 8.30676 11.13758',
        'C3 6 1.26895 9.02168 10.39032', 'C4 6 1.64225 10.07768 9.58845',
        'C5 6 2.98081 10.44432 9.51725', 'C6 6 3.92045 9.74974 10.25408',
        'C7 6 3.53891 8.69091 11.05301']
        >>> exp.format_atoms_for_export(cell, atoms, gui=True) # doctest: +NORMALIZE_WHITESPACE +REPORT_NDIFF
        ['C1 6 1.78099 7.14907 12.00423', 'C2 6 2.20089 8.30676 11.13758',
        'C3 6 1.26895 9.02168 10.39032', 'C4 6 1.64225 10.07768 9.58845',
        'C5 6 2.98081 10.44432 9.51725', 'C6 6 3.92045 9.74974 10.25408',
        'C7 6 3.53891 8.69091 11.05301']
        """
        el = Element()
        atoms = deepcopy(atoms)
        from misc import frac_to_cart
        for line in atoms:
            if int(line[1]) < 0:
                line[1] = int(line[1])  # abs(int(line[1])) <- No abs(), it causes confusion with GUI
            else:
                line[1] = el.get_atomic_number(el.get_atomlabel(line[0]))
            coord = frac_to_cart(line[2:5], cell)
            line[2:5] = coord
        newlist = []
        if not gui:
            for i in atoms:
                newlist.append('{:4.4s} {:4d} {:>7.4f}  {:>7.4f}  {:>7.4f}'.format(*i))
        else:
            for i in atoms:
                newlist.append('{} {} {:>7.5f} {:>7.5f} {:>7.5f}'.format(i[0], i[1], i[2], i[3], i[4]))
        return newlist

    def export_to_clip(self, fragname):
        fragname = fragname.lower()
        try:
            tst = self.copy_to_clipboard(fragname)
        except AttributeError as e:
            tst = False
            print(e)
        if tst:
            print('Exported "{0}" to the clipboard.'.format(fragname))
            return True
        else:
            return False

    def export_to_gui(self, fragname):
        """
        exports atoms to output for the DSRGui
        >>> gdb = ParseDB('../dsr_db.txt')
        >>> exp = Export(gdb=gdb)
        >>> print(exp.export_to_gui(fragname="toluene")) # doctest: +NORMALIZE_WHITESPACE +REPORT_NDIFF +ELLIPSIS
        C1 6 1.78099 7.14907 12.00423;;C2 6 2.20089 8.30676 11.13758;;C3 6 1.26895 9.02168 10.39032;;C4
        6 1.64225 10.07768 9.58845;;C5 6 2.98081 10.44432 9.51725;;C6 6 3.92045 9.74974 10.25408;;C7 6 3.53891 8.69091 11.05301
        >>> print(exp.export_to_gui(fragname="TOLUENE")) # doctest: +NORMALIZE_WHITESPACE +REPORT_NDIFF +ELLIPSIS
        C1 6 1.78099 7.14907 12.00423;;C2 6 2.20089 8.30676 11.13758;;C3 6 1.26895 9.02168 10.39032;;C4 6 1.64225 10.07768 9.58845;;C5 6 2.98081 10.44432 9.51725;;C6 6 3.92045 9.74974 10.25408;;C7 6 3.53891 8.69091 11.05301
        """
        fragname = fragname.lower()
        atoms = self.format_atoms_for_export(self._gdb.get_cell(fragname), self._gdb.get_atoms(fragname), gui=True)
        atoms = ';;'.join(atoms)
        return atoms

    def file_is_opened(self, base, ending):
        """
        determines if the filebase.ending is opened and locked
        returns True if file can be opened
        returns False if file is locked
        """
        if '.' not in ending:
            ending = '.' + ending
        arg = base + ending
        try:
            print('try to open...')
            open(arg, 'w')
            return True
        except IOError:
            print('can not open file')
            return False

    def write_res_file(self, fragment):
        """
        Writes the atom data to a "self._fragment_name".res file
        """
        fragment = fragment.lower()
        # write to file:
        resfile = str(fragment) + '.res'
        try:
            f = open(resfile, 'w')
            for line in self.export_resfile(fragment):
                f.write(''.join(line))
            print('Database entry of "{}" successfully written to {}.'.format(fragment, resfile))
        except IOError:
            print('*** Could not write file {} ***'.format(resfile))
            import sys
            sys.exit()
        f.close()


if __name__ == '__main__':

    import sys
    import doctest

    failed, attempted = doctest.testmod()  # verbose=True)
    if failed == 0:
        print('passed all {} tests!'.format(attempted))
    else:
        print('{} of {} tests failed'.format(failed, attempted))

    sys.exit()

    ##############################################################################

    #    reference.close()
    #    http://en.wikipedia.org/wiki/Molecular_graphics
    #    // assume:
    #    // atoms with x, y, z coordinates (Angstrom) and elementSymbol
    #    // bonds with pointers/references to atoms at ends (bond pairs)
    #    // need pairs of bonded atoms and thus covalence radii
    #    // table of colors for elementTypes
    #    // find limits of molecule in molecule coordinates as xMin, yMin, xMax, yMax
    #    scale = min(xScreenMax/(xMax-xMin), yScreenMax/(yMax-yMin))
    #    xOffset = -xMin * scale; yOffset = -yMin * scale
    #    for (bond in $bonds) {
    #        atom0 = bond.getAtom(0)
    #        atom1 = bond.getAtom(1)
    #        x0 = xOffset+atom0.getX()*scale
    #        y0 = yOffset+atom0.getY()*scale // (1)
    #        x1 = xOffset+atom1.getX()*scale
    #        y1 = yOffset+atom1.getY()*scale // (2)
    #        x1 = atom1.getX()
    #        y1 = atom1.getY()
    #        #xMid = (x0 + x1) /2;  yMid = (y0 + y1) /2;
    #        #color0 = ColorTable.getColor(atom0.getSymbol())
    #        drawLine (color0, x0, y0, xMid, yMid)
    #        #color1 = ColorTable.getColor(atom1.getSymbol())
    #        #drawLine (color1, x1, y1, xMid, yMid)
    #    }
