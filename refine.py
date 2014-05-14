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
import subprocess
import sys, os, re
import shutil
from resfile import ResList
import misc


__metaclass__ = type  # use new-style classes

class ShelxlRefine():
    '''A class to do a shelxl refinement. It is only for shelxl 2013!
    The resfilename should be without ending.
    '''
    def __init__(self, reslist, resfile_name, find_atoms):
        self._find_atoms = find_atoms
        self._atoms_in_reslist = self._find_atoms.collect_residues()
        self.atoms = []
        for residues in list(self._atoms_in_reslist.keys()):
            for i in self._atoms_in_reslist[residues]:
                self.atoms.append(i[0])
        self.number_of_atoms = len(self.atoms)
        self.resfile_name = str(resfile_name)
        self.__reslist = reslist
        self.__shelx_command = self.find_shelxl_exe()
        self.b_array = self.approx_natoms()
        self.bakfile = str(self.resfile_name+'.dsr-bak')
        if not self.__shelx_command:
            print('\nSHELXL executable not found! No fragment fitting possible.\n')
            print('You can download SHELXL at http://shelx.uni-ac.gwdg.de/SHELX/index.php\n')
            sys.exit()
    

    def get_xl_version_string(self, exe):
        with open(exe, 'rb') as f:
            binary = f.read()
            position = binary.find(b'Version 201')
            if position > 0:
                f.seek(position+8, 0) # seek to version string
                version = f.read(6)   # read version string
                return version.decode('ascii')


    def find_shelxl_exe(self):
        '''
        returns the appropriate shelxl executable
        filesize 2013/4: 4900352
        filesize 97    : 1703936
        '''
        names = ['shelxl', 'xl']
        shx_exe = []
        for name in names:
            shx_exe.extend(misc.which(name))  # list of shelxl executables in path
            try:
                exe = shx_exe[0]
            except(IndexError):
                continue
            version = self.get_xl_version_string(exe)
            if not version:
                print('Your SHELXL version', exe, 'is too old for this Program')
                print('Please use SHELXL 2013/4 or above!')
                sys.exit(-1)
            version = version.split('/')
            if int(version[0]) < 2013:
                print('Your SHELXL version is too old. Please use SHELXL 2013/4 or above!')
                print('You can download SHELXL at http://shelx.uni-ac.gwdg.de/SHELX')
                sys.exit()
            else:
                return exe
                
    
    def approx_natoms(self):
        '''
        approximates the number of atoms to set the B array accordingly
        it is slightly more because q-peaks are also counted as atoms
        '''
        barray = self.number_of_atoms * 9
        return barray
        
    
    def set_refinement_cycles(self, cycles='8'):
        '''
        Modifies the number of refinement cycles in the reslist. 
        '''
        status=self.checkFileExist(self.resfile_name+'.res')
        if not status:
            print('Error: unable to find res file!')
            sys.exit(-1)
        regex = '(^L\.S\.\s)|(^CGLS\s)'
        ls_line = misc.find_line(self.__reslist, regex)
        ls_list = self.__reslist[ls_line].split()
        ls_list[1] = str(cycles)+'\n'
        self.__reslist[ls_line] = '  '.join(ls_list)


        
    def remove_acta(self):
        '''
        removes the ACTA card, bacause we refine with L.S. 0, wich is incompatible!
        '''
        regex = '^ACTA'
        acta_line = misc.find_multi_lines(self.__reslist, regex)
        if acta_line:
            for i in acta_line:
                self.__reslist[i] = 'rem ACTA\n'
        
    
    def afix_is_closed(self, line):
        '''
        check if last afix was closed with afix 0
        
        - returns False if last AFIX was not closed
        - returns True if last AFIX was closed
        - returns True if no AFIX found at all
        '''
        afixes = []
        for num, i in enumerate(self.__reslist):
            i = i.upper()
            if num <= line:
                if i.startswith('AFIX'):
                    afixes.append(i.split()[1])
        try:
            if str(afixes[-2]) != '0':
                return False
            else:
                return True
        except(IndexError):
            return True


    
    def remove_afix(self):
        '''
        removes the AFIX 17X after refinement.
        note: find_line matches case insensitive
        '''
        regex = '^AFIX\s+9'
        afix_line = misc.find_line(self.__reslist, regex)
        if self.afix_is_closed(afix_line): # only delete afix if last afix was closed
            if afix_line:
                del self.__reslist[afix_line]
        else:
            self.__reslist[afix_line] = 'AFIX 0\n'
        
    

    def checkFileExist(self, filename):
        '''
        check if shelxl as written a res file successfully
        '''
        status = False
        res_filesize = 0
        try:
            res_filesize = int(os.stat(str(filename)).st_size)
        except:
            print('"{}" not found!'.format(filename))
            status = False
        if res_filesize > 10:
            status = True
        return status

        
    def backup_shx_file(self):
        '''
        makes a copy of the res file
        '''
        resfile = str(self.resfile_name+'.res')
        try:
            shutil.copyfile(resfile, self.bakfile)
        except(IOError):
            print('Unable to make backup file from {}.'.format(resfile))
            sys.exit(-1)


    def restore_shx_file(self):
        '''
        restores filename from backup
        '''
        resfile = str(self.resfile_name+'.res')
        try:
            shutil.copyfile(self.bakfile, resfile)
        except(IOError):
            print('Unable to make restore res file from {}.'.format(self.bakfile))
        try:
            misc.remove_file(self.bakfile)
        except(IOError):
            print('Unable to delete backup file {}.'.format(self.bakfile))


    def pretty_shx_output(self, out):
        '''
        selectively prints the output from shelx
        '''
        wr2 = False
        r1 = False
        gof = False
        for i in out:
            #if re.match(r'.*Command line parameters', i):
            if i.startswith(' +  Copyright(C)'):
                print(' SHELXL '+' '.join(i.split()[6:8]))
            # wR2
           # These values are always bad after a simple LS fit without any atom movement:
            if i.startswith(' wR2') and not wr2:
                wr2 = True
                line = i[:].split()
                print(' {}  {} {:>6}'.format(line[0], line[1], line[2][:6]))
            # R1
            if i.startswith(' R1') and not r1:
                r1 = True
                line = i[:].split()
                print(' {}   {} {:>6}'.format(line[0], line[1], line[2][:6]))
            # GooF
            if re.match(r'.*GooF.*', i) and not gof:
                gof = True
                line = i.split()
                print(' {} {} {:>5}0'.format(line[0], line[1], line[4][:5]))
            if re.match(r'.*CANNOT\s+OPEN\s+FILE.*hkl.*', i):
                print(' No hkl file found!')
                print('You need a proper hkl file to use DSR!')
                sys.exit()
            if re.match(r'.*\*\*.*', i):
                print(' SHELXL says:')
                print(' {}'.format(i.strip('\n\r')))
    

    def run_shelxl(self):
        '''
        This method runs shelxl 2013 on the res file self.resfile_name
        '''
        
        resfile = self.resfile_name+'.res'
        hklfile = self.resfile_name+'.hkl'
        
        if not self.checkFileExist(hklfile):
            print('You need a proper hkl file to use DSR.')
            sys.exit()
        
        command_line='{} -b{} {}'.format(self.__shelx_command, self.b_array, self.resfile_name).split()
        
        self.backup_shx_file()
        print('-----------------------------------------------------------------')
        print(' refining with "{}" and "L.S. 0"'.format(' '.join(command_line)))
        p = subprocess.Popen(command_line, stdin = subprocess.PIPE,
                            stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
        (child_stdin, child_stdout_and_stderr) = (p.stdin, p.stdout)
        child_stdin.close()
        
        # Watch the output for successful termination
        out = child_stdout_and_stderr.readline().decode('ascii')
        output = []
        while out:
            output.append(out)
            out = child_stdout_and_stderr.readline().decode('ascii')
        
        child_stdout_and_stderr.close()
        # output only the most importand things from shelxl:
        self.pretty_shx_output(output)

        status = self.checkFileExist(resfile) # status is False if shelx was unsecessful

        if not status: # fail
            print('-----------------------------------------------------------------')
            print('\nError: SHELXL terminated unexpectedly. Restoring original file.')
            print('Check for errors in your SHELX input file!\n')
            self.restore_shx_file()
            sys.exit()
        else:          # sucess
            print('-----------------------------------------------------------------')
            try:
                misc.remove_file(self.bakfile)
            except(IOError):
                print('Unable to delete backup file {}.'.format(self.bakfile))
            
        
        
        
if __name__ == '__main__':
    rl = ResList(res_file)
    res_list = rl.get_res_list()
    
    shx = ShelxlRefine(res_list, 'testfile')
    print(shx.afix_is_closed(110))
    print(shx.remove_afix())
    #shx.set_refinement_cycles()
    #shx.run_shelxl()
