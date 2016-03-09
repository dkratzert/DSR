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
from misc import find_line, checkFileExist, remove_line
from constants import sep_line



__metaclass__ = type  # use new-style classes

class ShelxlRefine():
    '''
    A class to do a shelxl refinement. It is only for shelxl 2013 and above!
    The resfilename should be without ending.
    '''
    def __init__(self, reslist, resfile_name, find_atoms):
        '''
        :param reslist: SHELXL .res file as list
        :param resfile_name: SHELXL res file name
        :param find_atoms: FindAtoms() object
        '''
        self._find_atoms = find_atoms
        self._atoms_in_reslist = self._find_atoms.collect_residues()
        self.number_of_atoms = 0
        for residues in self._atoms_in_reslist:
            self.number_of_atoms+=len(self._atoms_in_reslist[residues])
        self.resfile_name = str(resfile_name)
        self._reslist = reslist
        self._shelx_command = self.find_shelxl_exe()
        self.b_array = self.approx_natoms()
        self.backup_file = str(self.resfile_name+'.dsr-bak')

        if not self._shelx_command:
            print('\nSHELXL executable not found in system path! No fragment fitting possible.\n')
            print('You can download SHELXL at http://shelx.uni-ac.gwdg.de/SHELX/index.php\n')
            sys.exit()


    def get_xl_version_string(self, exe):
        """
        Extracts the version string from a SHELXL executable.
        This is fast and needs no hashes etc.
        """
        try:
            with open(exe, 'rb') as f:
                binary = f.read()
                position = binary.find(b'Version 201')
                if position > 0:
                    f.seek(position+8, 0) # seek to version string
                    version = f.read(6)   # read version string
                    return version.decode('ascii')
                else:
                    return None
        except(IOError):
            print("Could not determine SHELXL version. DSR might fail to run.")
            return None

    def find_shelxl_exe(self):
        '''
        returns the appropriate shelxl executable
        '''
        names = ['shelxl', 'xl']
        download = 'You can download SHELXL at http://shelx.uni-ac.gwdg.de/SHELX'
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
                print(download)
                sys.exit(-1)
            version = version.split('/')
            if int(version[0]) < 2013:
                print('Your SHELXL version is too old. Please use SHELXL 2013/4 or above!')
                print(download)
                sys.exit()
            else:
                return exe


    def approx_natoms(self):
        '''
        approximates the number of atoms to set the B array accordingly
        it is slightly more because q-peaks are also counted as atoms
        '''
        barray = self.number_of_atoms * 9
        if barray <= 3000:
            barray = 3000
        return barray

    def set_refinement_cycles(self, cycles):
        '''
        Modifies the number of refinement cycles in the reslist.
        '''
        status = checkFileExist(self.resfile_name + '.res')
        if not status:
            print('Error: unable to find res file!')
            sys.exit(-1)
        regex = r'(^L\.S\.\s)|(^CGLS\s)'
        ls_line = misc.find_line(self._reslist, regex)
        ls_list = self._reslist[ls_line].split()
        ls_list[1] = str(cycles)+'\n'
        self._reslist[ls_line] = '  '.join(ls_list)
        return self._reslist

    @property
    def get_refinement_cycles(self):
        '''
        Modifies the number of refinement cycles in the reslist.
        '''
        status = checkFileExist(self.resfile_name + '.res')
        if not status:
            print('Error: unable to find res file!')
            sys.exit(-1)
        regex = r'(^L\.S\.\s)|(^CGLS\s)'
        ls_line = misc.find_line(self._reslist, regex)
        try:
            cycles = self._reslist[ls_line].split()[1]
        except(IndexError):
            return None
        return cycles

    def remove_acta_card(self):
        '''
        removes the ACTA card, bacause we refine with L.S. 0, wich is incompatible!
        '''
        regex = r'^ACTA'
        acta_lines = misc.find_multi_lines(self._reslist, regex)
        commands = []
        if acta_lines:
            for i in acta_lines:
                commands.append(self._reslist[i])
                self._reslist[i] = 'REM '+self._reslist[i]
        return commands

    def restore_acta_card(self, lines):
        '''
        restores original acta cards.

        :param lines: list of lines with previous acta cards
        '''
        if not lines:
            return
        regex = r'^rem\s{1,5}ACTA'
        acta_lines = misc.find_multi_lines(self._reslist, regex)
        if not len(lines) == len(acta_lines):
            return
        if acta_lines:
            for n, i in enumerate(acta_lines):
                self._reslist[i] = lines[n]

    def afix_is_closed(self, line):
        '''
        check if last afix before dsr command was closed with afix 0

        - returns False if last AFIX was not closed
        - returns True if last AFIX was closed
        - returns True if no AFIX found at all
        '''
        afixes = []
        for num, i in enumerate(self._reslist):
            i = i.upper()
            if num <= line:
                if re.match(r'AFIX', i, re.IGNORECASE):
                    afixes.append(i.split()[1])
        try:
            if str(afixes[-2]) != '0':
                return False # last afix not closed
            else:
                return True # last afix is closed
        except(IndexError):
            # in this case, no other afix is present and deletion of afix 9 is save
            return True


    def remove_afix(self, random_num):
        '''
        removes the AFIX 9 after refinement.
        note: find_line matches case insensitive
        '''
        regex = 'REM '+random_num
        afix_line = misc.find_line(self._reslist, regex)
        if afix_line:
            #self._reslist[afix_line-1] = self._reslist[afix_line-1][4:]
            remove_line(self._reslist, afix_line, remove=True)
        # only delete afix if last afix before the dsr command was closed
        if self.afix_is_closed(afix_line):
            if afix_line:
                remove_line(self._reslist, afix_line-1, remove=True)
                afix_line2 = misc.find_line(self._reslist, regex)
                if afix_line2:
                    remove_line(self._reslist, afix_line2, remove=True)
                    remove_line(self._reslist, afix_line2, remove=True)
        else:
            if afix_line:
                self._reslist[afix_line-1] = 'AFIX 0\n'


    def backup_shx_file(self):
        '''
        makes a copy of the res file
        make backup in dsrsaves before every fragment fit.
        name: self.resfile_name-date-time-seconds.res
        '''
        import datetime
        bakup_dir = 'dsrsaves'
        now = datetime.datetime.now()
        timestamp = (str(now.year)+'_'+str(now.month)+'_'+str(now.day)+'_'+
                     str(now.hour)+'-'+str(now.minute)+'-'+str(now.second))
        resfile = str(self.resfile_name+'.res')
        try:
            shutil.copyfile(resfile, self.backup_file)
        except(IOError):
            print('Unable to make backup file from {}.'.format(resfile))
            sys.exit(-1)
        if not os.path.exists(bakup_dir):
            try:
                os.makedirs(bakup_dir)
            except(IOError, OSError):
                print('Unable to create backup directory {}.'.format(bakup_dir))
                #sys.exit(False)
        try:
            shutil.copyfile(resfile, bakup_dir+'/'+self.resfile_name+'_'+timestamp+'.res')
        except(IOError):
            print('\nUnable to make backup file from {} in dsrsaves.'.format(resfile))
            #sys.exit(-1)


    def restore_shx_file(self):
        '''
        restores filename from backup
        '''
        resfile = str(self.resfile_name+'.res')
        try:
            print('Restoring previous res file.')
            shutil.copyfile(self.backup_file, resfile)
        except(IOError):
            print('Unable to make restore res file from {}.'.format(self.backup_file))
        try:
            misc.remove_file(self.backup_file)
        except(IOError):
            print('Unable to delete backup file {}.'.format(self.backup_file))


    def pretty_shx_output(self, output):
        '''
        selectively prints the output from shelx
        '''
        wr2 = False
        r1 = False
        gof = False
        for out in output:
            #if re.match(r'.*Command line parameters', out):
            if out.startswith(' +  Copyright(C)'):
                print(' SHELXL '+' '.join(out.split()[6:8]))
            # wR2
            # These values are always bad after a simple LS fit without any atom movement:
            if out.startswith(' wR2') and not wr2:
                wr2 = True
                line = out[:].split()
                print(' {}  {} {:>6}'.format(line[0], line[1], line[2][:6]))
            # R1
            if out.startswith(' R1') and not r1:
                r1 = True
                line = out[:].split()
                print(' {}   {} {:>6}'.format(line[0], line[1], line[2][:6]))
            # GooF
            if re.match(r'.*GooF.*', out) and not gof:
                gof = True
                line = out.split()
                print(' {} {} {:>5}0'.format(line[0], line[1], line[4][:5]))
            if re.match(r'.*CANNOT RESOLVE (SAME|RIGU|SIMU|DELU)', out):
                print('\nWarning: Are you sure that all atoms are in the correct order?\n')
            if re.match(r'.*CANNOT\s+OPEN\s+FILE.*hkl.*', out):
                print(' No hkl file found!')
                print('You need a proper hkl file to use DSR!')
                sys.exit()    
            if re.match(r'.*\*\* Extinction \(EXTI\) or solvent.*', out):
                continue
            if re.match(r'.*\*\* MERG code changed to 0', out):
                # disable this output
                continue
            if re.match(r'.*\*\* Bond\(s\) to .* ignored', out):
                # disable this output
                continue
            if re.match(r'.*\*\*.*', out):
                print(' SHELXL says:')
                print(' {}'.format(out.strip('\n\r')))


    def run_shelxl(self):
        '''
        This method runs shelxl 2013 on the res file self.resfile_name
        '''
        resfile = self.resfile_name+'.res'
        hklfile = self.resfile_name+'.hkl'
        if not checkFileExist(hklfile):
            print('You need a proper hkl file to use DSR.')
            sys.exit()
        command_line='{} -b{} {}'.format(self._shelx_command, self.b_array, self.resfile_name).split()
        self.backup_shx_file()
        print(sep_line)
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
        status = checkFileExist(resfile) # status is False if shelx was unsecessful
        if not status: # fail
            print(sep_line)
            print('\nError: SHELXL terminated unexpectedly.')
            print('Check for errors in your SHELX input file!\n')
            self.restore_shx_file()
            sys.exit()



    def check_refinement_results(self, list_file):
        '''
        Does some checks if the refinement makes sense e.g. if the data to parameter
        ratio is in an acceptable range.
        :param list_file: SHELXL listing file
        :type list_file: list
        '''
        is_resfile_there = misc.checkFileExist(self.resfile_name+'.res')
        if is_resfile_there and is_resfile_there == 'zero':
            print('Something failed in SHELXL. Please check your .ins and .lst file!')
            self.restore_shx_file()
            try:
                misc.remove_file(self.backup_file)
            except(IOError):
                print('Unable to delete backup file {}.'.format(self.backup_file))
            sys.exit()
        if not is_resfile_there:
            print('Something failed in SHELXL. Please check your .ins and .lst file!')
            self.restore_shx_file()
            try:
                misc.remove_file(self.backup_file)
            except(IOError):
                print('Unable to delete backup file {}.'.format(self.backup_file))            
            sys.exit()
        regex_final = r' Final Structure Factor Calculation.*\n'
        final_results = find_line(list_file, regex_final)
        # find data and parameters:
        try:
            dataobj = re.search(r'[0-9]+\s+data', list_file[final_results+4])
            data = float(dataobj.group(0).split()[0])
            parameterobj = re.search(r'[0-9]+\s+parameters', list_file[final_results+4])
            parameters = float(parameterobj.group(0).split()[0])
        except(AttributeError):
            return False
        try:
            data_to_parameter_ratio = data/parameters
        except(ZeroDivisionError):
            return False
        lattline = find_line(list_file, r'^ LATT.*')
        centro = None
        if lattline:
            try:
                latt = int(list_file[lattline].split()[1])
            except(ValueError):
                latt = False
            if latt > 0:
                centro = True
            else:
                centro = False
        if centro == True and data_to_parameter_ratio < 10:
            print('Warning! The data / parameter ratio is getting low (ratio = {:.1f})! \
                    '.format(data_to_parameter_ratio))
        if centro == False and data_to_parameter_ratio < 7.5:
            print('Warning! The data / parameter ratio is getting low (ratio = {:.1f})! \
                    '.format(data_to_parameter_ratio))
        try:
            misc.remove_file(self.backup_file)
        except(IOError):
            print('Unable to delete backup file {}.'.format(self.backup_file))



if __name__ == '__main__':
    from atomhandling import FindAtoms
    from restraints import ListFile
    res_file = 'p21c.res'
    rl = ResList(res_file)
    res_list = rl.get_res_list()
    find_atoms = FindAtoms(res_list)
    #lf = ListFile('p21c')
    #lst_file = lf.read_lst_file()
    shx = ShelxlRefine(res_list, 'testfile', find_atoms)


    #print(shx.afix_is_closed(110))
    #print(shx.remove_afix())
    #shx.set_refinement_cycles()
    #shx.run_shelxl()
