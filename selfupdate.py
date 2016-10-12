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
    import urllib
    response = urllib.urlopen('http://www.xs3-data.uni-freiburg.de/data/version.txt')
    version = response.readline().strip()
    return version


def update_dsr():
    """
    Updates the running DSR to the current version on the web server.
    """
    version = get_current_dsr_version()
    from dsr import VERSION
    if int(VERSION) < int(version):
        print('*** Current available version of DSR is {}. Performing upate ***'.format(version))
        get_update_package(version)
        print('*** Finished updating ***')
        return True
    if int(VERSION) >= int(version):
        print('*** DSR is already up to date ***')
        return False


def overwrite_dir(root_src_dir, root_dst_dir, move=True):
    """
    Moves the content of scrdir over destdir and overwrites all files.

    :param src_dir: source directory
    :param dst_dir: target directory
    :return: True/False
    """
    for src_dir, dirs, files in os.walk(root_src_dir):
        dst_dir = src_dir.replace(root_src_dir, root_dst_dir, 1)
        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir)
        for file_ in files:
            src_file = os.path.join(src_dir, file_)
            dst_file = os.path.join(dst_dir, file_)
            if os.path.exists(dst_file):
                os.remove(dst_file)
            if move:
                shutil.move(src_file, dst_dir)
            else:
                shutil.copy2(src_file, dst_dir)
    return True


def get_update_package(version):
    """
    Downloads the current DSR distribution from the web server and
    updates the files.

    :type version: int or string

    Returns
    -------
    True/False
    """
    dsrdir = "/Applications/DSR" # TODO: add real DSR_DIR here
    import urllib
    response = urllib.urlopen('http://www.xs3-data.uni-freiburg.de/data/DSR-{}.tar.gz'.format(version))
    with tempfile.NamedTemporaryFile(delete=False) as tmpfile:
        tmpfile.write(response.read())
    tmpdir = tempfile.mkdtemp()  # a temporary directory
    try:
        with tarfile.open(tmpfile.name) as tarobj:
            tarobj.extractall(path=tmpdir)
    except tarfile.ReadError:
        print('*** Cound not get update from server. If this problem persists, please update manually! ***')
        return False
    os.remove(tmpfile.name)
    overwrite_dir(os.path.join(tmpdir, "DSR-{}".format(version)), dsrdir, move=False)
    shutil.rmtree(tmpdir, ignore_errors=True)  # cleanup the files
    return True





if __name__ == "__main__":
    #import sys
    #import doctest
    #failed, attempted = doctest.testmod()  # verbose=True)
    #if failed == 0:
    #    print('passed all {} tests!'.format(attempted))

    update_dsr()
