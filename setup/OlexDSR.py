import os
import sys
import shutil
import re
import subprocess
import olex
import olx
from olexFunctions import OlexFunctions
OV = OlexFunctions()

'''
To run this script, type spy.OlexDSR(help) in Olex2
This script does not work very well. Better use the makros in the "custom.xld"
'''

def dsr(command):
    inputfilename = OV.FileName()
    try:
        shutil.copyfile(inputfilename+'.ins', inputfilename+'.res')
        #os.popen(command)
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
        print p.stdout.read()
    except Exception, err:
        print "DSR failed to run"
        print "This is why: %s" %err
        return

def OlexDSR():
    print "Olex2 to DSR Linker"
    # initialisation step
    exe_name = "dsr.bat"
    if not exe_name:
        print 'The DSR executable could not be located, aborting'
        return
        
    #end of the initialisation
    inputfilename = OV.FileName()
    print "Input file is: ", inputfilename

    command = "dsr.bat -r {}.res".format(inputfilename)
    #dsr(command)
    try:
        dsr(command)
        shutil.copyfile(inputfilename+'.res', inputfilename+'.ins')
        OV.AtReap(inputfilename)
        return True
    except Exception, err:
        print "DSR gave up. This is why: %s" %err
        return

def reopen():
    return OV.FileName()
        
OV.registerFunction(OlexDSR)
OV.registerFunction(reopen)
