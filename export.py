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
import sys, re, os
from options import OptionsParser
import atomhandling as at
from misc import ll_to_string
from dbfile import global_DB
import pyperclip


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
    
    def __init__(self, options, fragment=False):
        self.__options = options
        self.__fragment = fragment
        if self.__options.export_fragment:
            self.__fragment = self.__options.export_fragment
        if self.__options.export_all:
            self.__fragment = fragment
        if not self.__fragment:
            print('Please run this object with command line parameter -e "fragment" !')
            sys.exit()
        self._gdb = global_DB()
        try:
            self.__db = self._gdb.build_db_dict()[self.__fragment.lower()]
        except(KeyError):
            print('Fragment "{}" was not found in the database!!'.format(self.__fragment))
            sys.exit()
        self._comment = self.__db['comment']
        self.__dbatoms = self.__db['atoms']
        self.__atomtypes = at.get_atomtypes(self.__dbatoms)
        self.__fragline = self.__db['fragline']
        self.__cell = self.__fragline[2:]
        self.__clipcell = self.__fragline[:]
        self.copy_to_clipboard()
        self.format_calced_coords()  # expands the cell of calculated structures
        self.__cell = '  '.join(self.__cell)
        self._comment_regex = '^REM .*$'.upper()
        print('Exporting "{0}" to {0}.res'.format(self.__fragment))
    
    
    
    def format_calced_coords(self):
        '''
        In calculated structure the cell is 1 1 1 90 90 90. Shelxle has problems with that when growing.
        So the cell is expanded to 10 10 10
        '''
        summe = None
        summe = int(sum(float(i) for i in self.__cell[0:3])) # this is to detect calculated structures
        if summe == 3:  # 1+1+1=3!
            for coord in range(2,5):  # x, y, z of coordinates
                for line in self.__dbatoms:   # for every atom line
                    num = float(line[coord])/10
                    line[coord] = "{:10.6f}".format(num)
            # now the new 10,10,10 cell:
            for n in range(0,3):
                self.__cell[n] = '10'
        
    
    def export_resfile(self):
        '''exports a .res file from a database entry to be viewed in a GUI'''
        sfac = []
        res_export = []
        for i in self.__atomtypes:       #build sfac table from atomtypes
            if i not in sfac:
                sfac.append(i)

        atlist = []
        for i in self.__atomtypes:       #atomtypes in the db_entry
            for y, x in enumerate(sfac):
                if x == i:
                    atlist.append(y+1)
    
        for n, i in enumerate(atlist):
            self.__dbatoms[n][1] = i
        
        # build the UNIT table:
        unit = []
        for i in sfac:
            unit.append('1 ')   #no matter what number

        ## Now put all infos together:
        for i in self.__dbatoms:
            i[0] = i[0]+' ' # more space for long atom names
            i.append('11.00   0.04')  #make it a full qualified atom line with occupancy and U value
            
        a = ll_to_string(self.__dbatoms)      #convert to string
        res_export.append('TITL '+self.__fragment+'\n')     #title card with fragment name
        res_export.append('CELL  0.71  '+self.__cell+'\n')   # the cell with wavelength
        res_export.append('LATT  -1\n')
        res_export.append('SFAC '+'  '.join(sfac)+'\n')
        res_export.append('UNIT '+' '.join(unit)+'\n')
        res_export.append('WGHT  0.1'+'\n')
        res_export.append('FVAR  1'+'\n')
        res_export.append(a)                         # the atoms
        res_export.append('\nHKLF 4\nEND\n')         # the end
        return res_export

            
    def copy_to_clipboard(self):
        '''
        copys the exported atoms to the clipboard including FRAG cel FEND commands
        '''
        from misc import frac_to_cart
        import copy
        clip_text = []
        cell = self.__cell[:]
        cell = [float(x) for x in cell]
        atoms = copy.deepcopy(self.__dbatoms)
        for line in atoms:
            coord = frac_to_cart(line[2:5], cell)
            line[2:5] = coord
        clip_text.append('FRAG')
        clip_text.append('\n'+ll_to_string(atoms))
        clip_text.append('\nFEND')
        text = ' '.join(clip_text)
        pyperclip.setcb(text)


    def file_is_opened(self, base, ending):
        '''
        determines if the filebase.ending is opened and locked
        returns True if file can be opened
        returns False if file is locked
        '''
        arg = base+ending
        try:
            print('try to open...')
            f = open(arg, 'w')
            return True
        except IOError:
            print('can not open file')
            return False
        

    def write_file(self):
        '''Writes the data to a res file'''
        ## write to file:
        resfile = str(self.__fragment)+'.res'
        try:
            f = open(resfile, 'w')
            for line in self.export_resfile(): 
                line = ''.join(line)
                f.write(line)
            print('Database entry of "{}" successfully written to {}.'.format(self.__fragment, resfile))
        except(IOError):
            print('could not write file {}'.format(resfile))
            sys.exit(-1)
        f.close()
        #if self.__options.export_all:
        self.make_image()

    
    def make_image(self, debug=False):
        import shutil
        import time
        import misc
        import subprocess
        '''
        ellipsoid plot of the molecule
        '''
        resfile = str(self.__fragment)+'.res'
        insfile = str(self.__fragment)+'.ins'
        info=None
        commandline = 'platon -O {}'.format(insfile).split()
        try:
            info = subprocess.STARTUPINFO()
            info.dwFlags = 1
            info.wShowWindow = 0
        except(AttributeError):
            pass
            
        misc.remove_file(insfile) #platon runs faster if no ins file is present!
        misc.remove_file(self.__fragment+'.png', exit_dsr=True)
        if misc.which('platon'):
            pass
        else:
            print('Could not write a .png image. No platon executable in path.')
            sys.exit()
        try:
            shutil.copyfile(resfile, insfile)
        except(IOError):
            print('unable to write .ins file for plotting!')
            sys.exit(-1)

        try:
            plat = subprocess.Popen(commandline, stdin = subprocess.PIPE,
                            stdout = subprocess.PIPE, stderr=subprocess.STDOUT, startupinfo=info)
            timeticks = 0
            psfile = self.__fragment+'.ps'
            while not os.path.isfile(psfile):
                timeticks = timeticks+1
                time.sleep(0.01)
                if timeticks > 1500:
                    print('Platon run took too long to execute.')
                    break
            size1 = os.stat(psfile).st_size
            size2 = 99999999
            timeticks = 0
            while size1 < size2: 
                timeticks = timeticks+1
                size2 = os.stat(psfile).st_size
                time.sleep(0.1)
                if timeticks > 30:
                    plat.terminate()
                    break
        except(WindowsError) as e:
            print('unable to run platon!', e)
            extensions = ('.bin', '.def', '.hkp', '.ins', '.pjn', '.lis', '.res', '.sar', '.sum', '.eld', '.out')
            plat.terminate()
            for i in extensions:
                misc.remove_file(self.__fragment+i)
            sys.exit()
        
        misc.remove_file('platon.out', terminate=plat)
        extensions = ('.lis', '.eld', '.def', '.pjn')
        for i in extensions:
            misc.remove_file(self.__fragment+i)
        misc.remove_file(insfile)
        
        # test for convert from ImageMagic
        if misc.which('montage'): # i check for montage, because windows also ha a convert.exe
            pass
        else:
            print('Could not write a .png image. ImageMagic is not installed.')
            plat.terminate()
            return
        try:
            convert = 'convert'
            options_convert = '-crop 84%x90%+40%+40% -rotate 90 -trim' 
            print('converting from .ps to .png')
            files = '"{}.ps" "{}.png"'.format(self.__fragment, self.__fragment)
            image_commandline = '{} {} {}'.format(convert, options_convert, files)
            conv = os.popen(image_commandline)
            conv.close()
            # were we successful?
            if os.path.isfile(self.__fragment+'.png'):
                print('success!')
                # in case of success remove the postscript file
                misc.remove_file(psfile, terminate=plat)
            else:
                print('unable to write .png file. Is platon and ImageMagic installed?')
                plat.terminate()
        except(EnvironmentError) as e:
            print('unable to convert file', e)
        
        misc.remove_file(self.__fragment+'.lis')
        misc.remove_file(self.__fragment+'.eld')
        plat.terminate()
    
    
    



if __name__ == '__main__':
    from dbfile import global_DB
    gdb = global_DB()
    db = gdb.build_db_dict()['toluene']
    
    export = Export(options, 'toluene')
    for i in export.export_resfile():
        print(i.strip('\n'))
    #import pyperclip
    #pyperclip.setcb('The text to be copied to the clipboard.')
    #spam = pyperclip.getcb()
    #export.write_file()
    #export.make_image(debug=True)
    
