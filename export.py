# -*- encoding: utf-8 -*-
# möp
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

import copy
import os

import atomhandling as at
from atoms import Element
from dbfile import global_DB
from misc import wrap_headlines
from restraints import Restraints

__metaclass__ = type  # use new-style classes


class Export():
    '''
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
    '''

    def __init__(self, gdb, invert=False):
        '''

        :param fragment_name: string, name of the database fragment
        :param invert:        bool, should the coordinates be inverted?
        '''
        self.invert = invert
        self._gdb = gdb

    def format_calced_coords(self, cell, fragment):
        """
        TODO: make this crude hack more elegant!
        In calculated structure the cell is 1 1 1 90 90 90. Shelxle has problems
        with that when growing. So the cell is expanded to 50 50 50
        >>> gdb = global_DB()
        >>> exp = Export(gdb, invert=False)
        >>> exp.format_calced_coords([1, 1, 1, 90, 90, 90], "benzene")
        [['50', '50', '50', 90, 90, 90], [['C1', '1', '  0.017600', ' -0.006618', '  0.005344'], ['C2', '1', '  0.015724', ' -0.007554', '  0.004762'], ['C3', '1', '  0.015212', ' -0.006368', '  0.003851'], ['C4', '1', '  0.016584', ' -0.004244', '  0.003510'], ['C5', '1', '  0.018464', ' -0.003288', '  0.004080'], ['C6', '1', '  0.018976', ' -0.004464', '  0.004998']]]
        >>> exp.format_calced_coords([1, 1, 1, 90, 90, 90], "BENZENE")
        [['50', '50', '50', 90, 90, 90], [['C1', '1', '  0.017600', ' -0.006618', '  0.005344'], ['C2', '1', '  0.015724', ' -0.007554', '  0.004762'], ['C3', '1', '  0.015212', ' -0.006368', '  0.003851'], ['C4', '1', '  0.016584', ' -0.004244', '  0.003510'], ['C5', '1', '  0.018464', ' -0.003288', '  0.004080'], ['C6', '1', '  0.018976', ' -0.004464', '  0.004998']]]
        """
        fragment = fragment.lower()
        atoms = copy.deepcopy(self._gdb[fragment]['atoms'])
        summe = int(sum(float(i) for i in cell[0:3]))  # this is to detect calculated structures
        if summe == 3:  # 1+1+1=3!
            for coord in range(2, 5):  # x, y, z of coordinates
                for line in atoms:  # for every atom line
                    num = float(line[coord]) / 50
                    line[coord] = "{:10.6f}".format(num)
            # now the new 50,50,50 cell:
            for n in range(0, 3):
                cell[n] = '50'
        return [cell, atoms]

    def make_dfix(self, fragname):
        fragname = fragname.lower()
        restr = Restraints(self, fragname, self._gdb)
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
        for fragment in self._gdb.db_dict:
            self.write_res_file(fragment)
        sys.exit(1)

    def export_resfile(self, fragname):
        """
        exports a .res file from a database entry to be viewed in a GUI
        #>>> invert = False
        #>>> gdb = global_DB(invert)
        #>>> fragment = 'toLuene'
        #>>> export = Export(gdb, invert)
        """
        fragname = fragname.lower()
        comment = self._gdb[fragname]['comment']
        cell = self._gdb[fragname]['fragline'][2:]
        # expands the cell of calculated structures:
        cell, atoms = self.format_calced_coords(cell=cell, fragment=fragname)
        cellstring = ' {:>8.3f} {:>8.3f} {:>8.3f} {:>8.3f} {:>8.3f} {:>8.3f}'.format(*[float(i) for i in cell])
        print('Exporting "{0}" to {0}.res'.format(fragname))
        if self.invert:
            print("Fragment inverted.")
        try:
            from dsr import VERSION
        except(ImportError):
            VERSION = ''
        sfac = []
        res_export = []
        self._gdb.check_consistency(fragname)
        self._gdb.check_db_atom_consistency(fragname)
        for i in at.get_atomtypes(atoms):  # build sfac table from atomtypes
            if i not in sfac:
                sfac.append(i)

        atlist = []
        for i in at.get_atomtypes(atoms):  # atomtypes in the db_entry
            for y, x in enumerate(sfac):
                if x == i:
                    atlist.append(y + 1)

        for n, i in enumerate(atlist):
            atoms[n][1] = i

        # build the UNIT table:
        unit = []
        for i in sfac:
            unit.append('1 ')  # no matter what number

        ## Now put all infos together:
        for i in atoms:
            i[0] = i[0] + ' '  # more space for long atom names
            i.append('11.00   0.04')  # make it a full qualified atom line with occupancy and U value

        final_atomlist = [('{:4.4s} {:4.2s} {:>8.5f}  {:>8.5f}  {:>8.5f}   11.0   0.04\n'.format(
            str(i[0]), str(i[1]), float(i[2]), float(i[3]), float(i[4]))) for i in atoms]
        res_export.append('TITL ' + fragname + '\n')  # title card with fragment name
        try:
            res_export.append('REM This file was exported by DSR version {}\n'.format(VERSION))
        except(NameError):
            pass
        res_export.append('REM ' + '\nREM '.join(comment) + '\n')
        res_export.append('CELL 0.71073 ' + cellstring + '\n')  # the cell with wavelength
        res_export.append('ZERR    1.00   0.000    0.000    0.000    0.000    0.000    0.000\n')
        res_export.append('LATT  -1\n')
        res_export.append('SFAC ' + '  '.join(sfac) + '\n')
        res_export.append('UNIT ' + ' '.join(unit) + '\n')
        res_export.append('REM  RESIDUE: {}\n'.format(self._gdb.get_resi_from_fragment(fragname)))
        res_export.append('REM Sum formula: {}\n'.format(self._gdb.get_sum_formula(fragname)))
        res_export.append('WGHT  0.1' + '\n')
        res_export.append('FVAR  1.0' + '\n')
        try:
            res_export.append('rem Restraints from DSR database:\n')
            res_export.append(''.join(wrap_headlines(self._gdb.get_head_from_fragment(fragname))))
        except:
            pass
        try:
            res_export.append('rem Restraints from atom connectivities:\n')
            res_export.append(self.make_dfix(fragname))
            res_export.append('rem end of restraints\n')
        except:
            pass
        res_export.append('\n')
        res_export.append(final_atomlist)  # the atoms
        res_export.append('\nHKLF 0\nEND\n')  # the end
        return res_export

    def copy_to_clipboard(self, fragname):
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
        atoms = self.format_atoms_for_export(fragname)
        atoms = '\n'.join(atoms)
        clip_text.append('FRAG')
        clip_text.append('\n' + atoms)
        clip_text.append('\nFEND')
        text = ' '.join(clip_text)
        pyperclip.setcb(text)
        return True

    def format_atoms_for_export(self, fragname, gui=False):
        """
        fractional coordinates are converted to cartesian
        Atom;;number;;x;;y;;z

        >>> gdb = global_DB(invert=False)
        >>> exp = Export(gdb=gdb, invert=False)
        >>> print(exp.format_atoms_for_export(fragname="toluene", gui=False)) # doctest: +NORMALIZE_WHITESPACE +REPORT_NDIFF
        ['C1      6  1.7810   7.1491  12.0042',
        'C2      6  2.2009   8.3068  11.1376',
        'C3      6  1.2689   9.0217  10.3903',
        'C4      6  1.6422  10.0777   9.5884',
        'C5      6  2.9808  10.4443   9.5173',
        'C6      6  3.9205   9.7497  10.2541',
        'C7      6  3.5389   8.6909  11.0530']

        >>> print(exp.format_atoms_for_export(fragname="toluene", gui=True)) # doctest: +NORMALIZE_WHITESPACE +REPORT_NDIFF
        ['C1 6 1.78099 7.14907 12.00423', 'C2 6 2.20089 8.30676 11.13758',
        'C3 6 1.26895 9.02168 10.39032', 'C4 6 1.64225 10.07768 9.58845',
        'C5 6 2.98081 10.44432 9.51725', 'C6 6 3.92045 9.74974 10.25408',
        'C7 6 3.53891 8.69091 11.05301']
        
        >>> print(exp.format_atoms_for_export(fragname="TOLUeNE", gui=True)) # doctest: +NORMALIZE_WHITESPACE +REPORT_NDIFF
        ['C1 6 1.78099 7.14907 12.00423', 'C2 6 2.20089 8.30676 11.13758',
        'C3 6 1.26895 9.02168 10.39032', 'C4 6 1.64225 10.07768 9.58845',
        'C5 6 2.98081 10.44432 9.51725', 'C6 6 3.92045 9.74974 10.25408',
        'C7 6 3.53891 8.69091 11.05301']
        
        #>>> print(exp.format_atoms_for_export(fragname="", gui=True))
        #sys.exit()
        """
        fragname = fragname.lower()
        el = Element()
        from misc import frac_to_cart
        dbentry = self._gdb[fragname]
        cell = dbentry['fragline'][2:]
        cell = [float(x) for x in cell]
        atoms = copy.deepcopy(dbentry['atoms'])
        for line in atoms:
            if int(line[1]) < 0:
                line[1] = int(line[1])  # abs(int(line[1])) <- No abs(), it causes confusion with GUI
            else:
                line[1] = el.get_atomic_number(el.get_atomlabel(line[0]))
            frac_coord = [float(i) for i in line[2:5]]
            coord = frac_to_cart(frac_coord, cell)
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
        except(AttributeError) as e:
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

        >>> gdb = global_DB(invert=False)
        >>> exp = Export(gdb=gdb, invert=False)
        >>> print(exp.export_to_gui(fragname="toluene")) # doctest: +NORMALIZE_WHITESPACE +REPORT_NDIFF +ELLIPSIS
        C1 6 1.78099 7.14907 12.00423;;C2 6 2.20089 8.30676 11.13758;;C3 6 1.26895 9.02168 10.39032;;C4
        6 1.64225 10.07768 9.58845;;C5 6 2.98081 10.44432 9.51725;;C6 6 3.92045 9.74974 10.25408;;C7 6 3.53891 8.69091 11.05301
        >>> print(exp.export_to_gui(fragname="TOLUENE")) # doctest: +NORMALIZE_WHITESPACE +REPORT_NDIFF +ELLIPSIS
        C1 6 1.78099 7.14907 12.00423;;C2 6 2.20089 8.30676 11.13758;;C3 6 1.26895 9.02168 10.39032;;C4 6 1.64225 10.07768 9.58845;;C5 6 2.98081 10.44432 9.51725;;C6 6 3.92045 9.74974 10.25408;;C7 6 3.53891 8.69091 11.05301
        """
        fragname = fragname.lower()
        atoms = self.format_atoms_for_export(fragname, gui=True)
        atoms = ';;'.join(atoms)
        return atoms

    def file_is_opened(self, base, ending):
        """
        determines if the filebase.ending is opened and locked
        returns True if file can be opened
        returns False if file is locked
        """
        if not '.' in ending:
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
                line = ''.join(line)
                f.write(line)
            print('Database entry of "{}" successfully written to {}.'.format(fragment, resfile))
        except IOError:
            print('*** Could not write file {} ***'.format(resfile))
            import sys
            sys.exit(-1)
        f.close()
        print("Image creation is currently disabled, because PLATON and ImageMagic are causing troubles.")
        #self.make_image(fragment)

    def make_image(self, fragname):
        from shutil import copyfile
        import time
        import misc
        import subprocess
        '''
        Draws an ellipsoid plot of the molecule. This method depends on PLATON 
        from Ton Spek. The windows version of Platon needs some special care 
        because of its nasty output window.
        This method tries to kill the platon process and removes all leftorver 
        file in case something goes wrong during the image drawing process.
        '''
        fragname = fragname.lower()
        plat = None
        resfile = str(fragname) + '.res'
        insfile = str(fragname) + '.ins'
        info = None
        commandline = 'platon -O {}'.format(insfile).split()
        try:
            info = subprocess.STARTUPINFO()
            info.dwFlags = 1
            info.wShowWindow = 0
        except(AttributeError):
            pass
        misc.remove_file(insfile)  # platon runs faster if no ins file is present!
        misc.remove_file(fragname + '.png', exit_dsr=True)
        if not misc.which('platon'):
            print('*** Could not write a .png image. No PLATON executable in PATH found. ***')
            return None
        try:
            copyfile(resfile, insfile)
        except(IOError):
            print('*** Unable to write .ins file for plotting! ***')
            return None
        try:
            plat = subprocess.Popen(commandline, stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, startupinfo=info)
            timeticks = 0
            psfile = fragname + '.ps'
            while not os.path.isfile(psfile):
                timeticks = timeticks + 1
                time.sleep(0.01)
                # give PLATON 15s to draw the picture
                if timeticks > 1500:
                    print('PLATON run took too long to execute. Killing Platon...')
                    try:
                        plat.terminate()
                    except:
                        sys.exit()
                    break
            size1 = os.stat(psfile).st_size
            size2 = 99999999
            timeticks = 0
            while size1 < size2:
                timeticks = timeticks + 1
                size2 = os.stat(psfile).st_size
                time.sleep(0.1)
                # give the system 3s to store the picture
                if timeticks > 30:
                    try:
                        plat.terminate()
                    except:
                        pass
                    break
        except() as e:
            print('unable to run platon!', e)
            extensions = ('.bin', '.def', '.hkp', '.ins', '.pjn', '_pl.spf',
                          '.lis', '.res', '.sar', '.sum', '.eld', '.out')
            try:
                plat.terminate()
            except:
                pass
            for i in extensions:
                # clean all the leftover files
                misc.remove_file(fragname + i)
            sys.exit()
        misc.remove_file('platon.out', terminate=plat)
        extensions = ('.lis', '.eld', '.def', '.pjn', '_pl.spf')
        for i in extensions:
            misc.remove_file(fragname + i)
        misc.remove_file(insfile)
        # test for convert from ImageMagic
        plat.terminate()


    def convert_eps_to_png(self, plat, fragname, psfile):
        """
        Convert to png
        ImageMagic from APEX causes problems
        """
        fragname = fragname.lower()
        import misc
        if misc.which('montage'):  # i check for montage, because windows also ha a convert.exe
            pass
        else:
            print('Could not write a .ps and .png image. ImageMagic is not installed.')
            plat.terminate()
            return
        try:
            convert = 'convert'
            options_convert = '-crop 84%x90%+40%+40% -rotate 90 -trim'
            print('converting from .ps to .png')
            files = '"{}.ps" "{}.png"'.format(fragname, fragname)
            image_commandline = '{} {} {}'.format(convert, options_convert, files)
            conv = os.popen(image_commandline)
            conv.close()
            # were we successful?
            if os.path.isfile(fragname + '.png'):
                print('success!')
                # in case of success remove the postscript file
                misc.remove_file(psfile, terminate=plat)
            else:
                print('Unable to write .png file. Is PLATON and ImageMagic installed?')
                plat.terminate()
                misc.remove_file(fragname + '.lis')
                misc.remove_file(fragname + '.eld')
                misc.remove_file(fragname + '_pl.spf')
        except(EnvironmentError) as e:
            print('unable to convert postscript file', e)
        misc.remove_file(fragname + '.lis')
        misc.remove_file(fragname + '.eld')
        misc.remove_file(fragname + '_pl.spf')


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




    # from dbfile import global_DB
    gdb = global_DB()
    db = gdb.build_db_dict()['toluene']

    # export = Export('toluene')
    # export.export_to_clip()

    from pngcanvas import PNGCanvas

    BUFSIZE = 8 * 1024  # Taken from filecmp module
    HEIGHT = WIDTH = 512
    c = PNGCanvas(WIDTH, HEIGHT, color=(0xff, 0, 0, 0xff))
    c.rectangle(0, 0, WIDTH - 1, HEIGHT - 2)
    c.rectangle(100, 100, 10, 10)
    c.filled_rectangle(100, 100, 10, 10)
    c.color = bytearray((0, 0, 0, 0xff))
    c.line(0, 0, WIDTH - 1, HEIGHT - 1)
    c.line(50, 50, 50, HEIGHT - 30)
    #      .|------|
    # atom1--|------|atom2
    #       |------|
    with open('reference.png', 'wb+') as reference:
        reference.write(c.dump())

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




    # for i in export.export_resfile():
    #    print(i.strip('\n'))
    # import pyperclip
    # pyperclip.setcb('The text to be copied to the clipboard.')
    # spam = pyperclip.getcb()
    # export.write_res_file()
    # export.make_image(debug=True)
