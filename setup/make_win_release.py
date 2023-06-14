#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This script has to be run from the main dir e.g. D:\GitHub\StructureFinder
"""
import pathlib
import subprocess
import sys

pth = pathlib.Path(__file__).parent.parent
print(pth)
sys.path.insert(0, str(pth / 'src/DSR'))
sys.path.insert(0, str(pth))

from version import VERSION

print("Updating version numbers to version {} ...".format(VERSION))

#print("Linux deb... {}".format(VERSION))
#process_debian_and_spec(debianpath)

#print("windows iss... {}".format(VERSION))
#process_iss(isspath)

print("Version numbers updated.")


def make_distribs():
    innosetup_compiler = r'C:\Program Files (x86)\Inno Setup 6/ISCC.exe'
    # Run DSR setup compiler
    subprocess.call([innosetup_compiler, f'/dMyAppVersion={VERSION}', r'setup\dsr-install.iss', ])


# Make a zip file for web interface distribution:
# make_zip(files)

# Make binary distributions:
make_distribs()
