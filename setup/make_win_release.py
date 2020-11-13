#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This script has to be run from the main dir e.g. D:\GitHub\StructureFinder
"""
import os
import subprocess
import sys

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from dsr import VERSION
from scripts.make_zipfile import make_zip, files
from setup.version_numbers import process_debian_and_spec, process_iss, debianpath, specpath, isspath

print("Updating version numbers to version {} ...".format(VERSION))

print("Linux deb... {}".format(VERSION))
process_debian_and_spec(debianpath)

print("Linux rpm... {}".format(VERSION))
process_debian_and_spec(specpath)

print("windows iss... {}".format(VERSION))
process_iss(isspath)

print("Version numbers updated.")


def make_distribs():
    innosetup_compiler = r'C:\Program Files (x86)\Inno Setup 6/ISCC.exe'
    # Run DSR setup compiler
    subprocess.call([innosetup_compiler, r'D:\GitHub\DSR\setup\dsr-install.iss', ])


# Make a zip file for web interface distribution:
make_zip(files)

# Make binary distributions:
make_distribs()

# This is not needed anymore:
# This is to test if the distribution numpy works:
#subprocess.call([r'C:\Python27-dsr\python.exe', '-c', 'import numpy as np; np.sin(3)'])

