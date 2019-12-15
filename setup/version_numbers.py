# /usr/bin/env python
# -*- encoding: utf-8 -*-

from __future__ import print_function

from dsr import VERSION
import misc
import dbfile
import os

"""
This file is for updating the various version number definitions for each DSR release
"""

debianpath = os.path.abspath("./setup/debian-package/DEBIAN/control")  # for the .deb file
isspath = os.path.abspath("./setup/dsr-install.iss")  # for the windows executable
specpath = os.path.abspath("./setup/dsr-linux.spec")  # for the .rpm file


def process_debian_and_spec(filepath):
    """
    search for: Version: 203
    replace with "Version: VERSION"
    """
    deb_file = dbfile.read_file_data(filepath)
    for num, line in enumerate(deb_file):
        if line.startswith("Version:"):
            l = line.split()
            l[1] = VERSION
            deb_file[num] = " ".join(l)
            break
    deb_file = [x + '\n' for x in deb_file]
    misc.write_file(deb_file, filepath)


def process_iss(filepath):
    """
    #define MyAppVersion  "203"
    """
    iss_file = dbfile.read_file_data(filepath)
    for num, line in enumerate(iss_file):
        if line.startswith("#define MyAppVersion"):
            l = line.split()
            l[2] = '"{}"'.format(VERSION)
            iss_file[num] = " ".join(l)
            break
    iss_file = [x + '\n' for x in iss_file]
    misc.write_file(iss_file, filepath)


if __name__ == "__main__":
    print("Updating version numbers to version {} ...".format(VERSION))

    print("Linux... {}".format(VERSION))
    process_debian_and_spec(debianpath)

    print("Linux... {}".format(VERSION))
    process_debian_and_spec(specpath)

    print("windows... {}".format(VERSION))
    process_iss(isspath)

    print("Version numbers updated.")
