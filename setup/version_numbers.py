#/usr/bin/env python
#-*- encoding: utf-8 -*-

from __future__ import print_function

from dsr import VERSION
import misc
import dbfile
import os

"""
This file is for updating the various version number definitions for each DSR release
"""

debianpath = os.path.abspath("./debian-package/DEBIAN/control")  # for the .deb file
isspath = os.path.abspath("./dsr-install.iss")  # for the windows executable
specpath = os.path.abspath("./dsr-linux.spec")  # for the .rpm file


def process_debian_and_spec(filepath):
    """
    search for: Version: 203
    replace with "Version: VERSION"
    """
    debfile = dbfile.read_file_data(filepath, with_comments=True)
    for num, line in enumerate(debfile):
        if line.startswith("Version:"):
            l = line.split()
            l[1] = VERSION
            debfile[num] = " ".join(l)+"\n"
            break
    misc.write_file(debfile, filepath)


def process_iss(filepath):
    """
    #define MyAppVersion  "203"
    """
    debfile = dbfile.read_file_data(filepath, with_comments=True)
    for num, line in enumerate(debfile):
        if line.startswith("#define MyAppVersion"):
            l = line.split()
            l[2] = '"{}"'.format(VERSION)
            debfile[num] = " ".join(l)+"\n"
            break
    misc.write_file(debfile, filepath)


if __name__ == "__main__":

    print("Updating version numbers...")

    print("Linux...")
    process_debian_and_spec(debianpath)

    print("Linux...")
    process_debian_and_spec(specpath)

    print("windows...")
    process_iss(isspath)

    print("Version numbers updated.")

