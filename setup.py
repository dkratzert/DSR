#/usr/bin/env python
#-*- encoding: utf-8 -*-
#möp

from distutils.core import setup

#import sys,os
"""
origIsSystemDLL = py2exe.build_exe.isSystemDLL()
def isSystemDLL(pathname):
        if os.path.basename(pathname).lower() in ("USER32.dll", "SHELL32.dll" ,"ADVAPI32.dll", "WS2_32.dll", "GDI32.dll", "KERNEL32.dll"):
                return 0
        return origIsSystemDLL(pathname)
py2exe.build_exe.isSystemDLL = isSystemDLL
"""



#This is a list of files to install, and where
#(relative to the 'root' dir, where setup.py is)
#You could be more specific.
files = ["*"]

setup(name = "DSR",
    version = "1.5.13",
    description = "yadda yadda",
    author = "Daniel Kratzert",
    author_email = "daniel.kratzert@ac.uni-freiburg.de",
    url = "https://www.xs3.uni-freiburg.de/research/dsr",
    packages = ['.'],
    #'package' package must contain files (see list above)
    #I called the package 'package' thus cleverly confusing the whole issue...
    #This dict maps the package name =to= directories
    #It says, package needs these files.
    package_data = {'package':  files },
    #'runner' is in the root.
    #scripts = ["runner"],
    long_description = """The program consists of a text database with fragments of molecules and the DSR program. 
    It acts as a preprocessor for SHELXL .res files. The user inserts a special command in the SHELXL .res file 
    and the DSR program reads this information to put a molecule or fragment with the desired atoms on the position 
    of the target atoms or q-peaks in the unit cell. Bond restraints can be either applied from the database to the molecule 
    or automatically generated.""" 
    #
    #This next part it for the Cheese Shop, look a little down the page.
    #classifiers = []     
) 