from sys import version_info
import sys
import os, stat
import re
from os.path import expanduser
import subprocess
import getpass

#if version_info > (2,6):
#    print 'DSR only works with Python version 2.6 and up. And also not python 3.x'
#    print 'Please install Python >2.6 first'
#    sys.exit()


home = expanduser("~")
dsr_dir='/opt/DSR'
dsr_path = 'PATH=$PATH:'+dsr_dir+'\nexport PATH'
DSR_DB_DIR = 'DSR_DB_DIR='+dsr_dir+'\nexport DSR_DB_DIR'
PROFILE = home+'/.profile'


try:
    os.access(dsr_dir, os.W_OK)
except:
    print 'no write permission'

def already_set_variable(file):
    for line in file:
        if 'opt/DSR' in line:
            return True
        else:
            return False
try:
#    os.chmod(dsr_dir, stat.S_IRWXU)
#    os.chmod(dsr_dir, stat.S_IRWXG)
    os.chmod(dsr_dir+"/dsr", stat.S_IRWXU)
    os.chmod(dsr_dir+"/dsr", stat.S_IRWXG)
except(OSError):
    print 'Please copy the content of the DSR dsr_linux-version.tar.gz file to '+dsr_dir+' first!'
    sys.exit()

with open(PROFILE, 'a+') as file:
    if not already_set_variable(file):
        print 'Environment path variable successfully set.'
        file.write(dsr_path+'\n')
        file.write(DSR_DB_DIR+'\n')
    else:
        print 'Environment path variable was already set.'
        sys.exit()


