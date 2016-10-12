# -*- encoding: utf-8 -*-
# möp
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

import shutil
import tarfile
import tempfile
import os

import misc


def get_current_dsr_version():
    """
    determines the current version of DSR on the web server

    >>> get_current_dsr_version()
    123

    Returns
    -------
    version number
    :type: int
    """
    pass

def update_dsr():
    """
    Updates the running DSR to the current version on the web server.

    >>> update_dsr()
    True

    Returns
    -------
    Sucess or not
    True/False
    """
    pass

def is_update_needed():
    """
    Decides if an update of DSR is needed and updates if neccesary. It does
    nothing in case of an already updated version.

    >>> is_update_needed()
    "updating DSR ..."
    "update was sucessful"
    "You are already using the current version. No update needed."

    Returns
    -------
    True/False
    """
    pass


def get_update_package(version, name):
    """
    Downloads the current DSR distribution from the web server and
    returns True if it suceeded.

    :type version: int or string
    :type name: string

    Returns
    -------
    True/False
    """
    import urllib
    response = urllib.urlopen('http://www.xs3-data.uni-freiburg.de/data/DSR-192.tar.gz')
    with tempfile.NamedTemporaryFile(delete=False) as tmpfile:
        tmpfile.write(response.read())
    tmpdir = tempfile.mkdtemp()
    with tarfile.open(tmpfile.name) as tarobj:
        tarobj.extractall(path=tmpdir)
    misc.remove_file(tmpfile.name)
    shutil.move(os.path.join(tmpdir, "DSR-192"), os.path.join(tmpdir, "DSR"))
    shutil.copy(os.path.join(tmpdir, "DSR"), "D:/Programme/")
    #shutil.rmtree(os.path.join(tmpdir, "DSR"))
    tmpdir.clear()




if __name__ == "__main__":
    import sys
    import doctest
    #failed, attempted = doctest.testmod()  # verbose=True)
    #if failed == 0:
    #    print('passed all {} tests!'.format(attempted))

    get_update_package('er', 'drsr')