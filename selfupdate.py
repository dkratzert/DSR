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

import os
import shutil
import tarfile
import tempfile
import urllib
import sys

from dsr import VERSION

urlprefix = "http://www.xs3-data.uni-freiburg.de/data"

# changes the user-agent of the http request:
class DSRURLopener(urllib.FancyURLopener):
    version = "DSR-updater"

urllib._urlopener = DSRURLopener()


def get_current_dsr_version(silent=False):
    """
    determines the current version of DSR on the web server

    >>> get_current_dsr_version()
    '195'

    Returns
    -------
    version number
    :type: int
    """
    try:
        response = urllib.urlopen('{}/version.txt'.format(urlprefix))
    except IOError:
        if not silent:
            print("*** Unable to connect to update server. No Update possible. ***")
        return 0
    version = response.readline().strip()
    return version


def is_update_needed(silent=False):
    """
    Decides if an update of DSR is needed.
    :return: True/False
    >>> is_update_needed()
    False
    """
    version = get_current_dsr_version(silent)
    if int(VERSION) < int(version):
        return True
    else:
        return False


def update_dsr(force=False, version=None):
    """
    Updates the running DSR to the current version on the web server.
    """
    if version:
        version = version
    else:
        version = get_current_dsr_version()
    if force:
        status = get_update_package(version)
        if status:
            print('*** Finished updating to version {} ***'.format(version))
            return True
        else:
            return False
    if int(VERSION) < int(version):
        print('*** Current available version of DSR is {}. Performing upate ***'.format(version))
        status = get_update_package(version)
        if status:
            print('*** Finished updating to version {} ***'.format(version))
            return True
        else:
            return False
    if int(VERSION) >= int(version):
        print('*** DSR is already up to date (version {}) ***'.format(version))
        return False


def get_system():
    """
    Gets the currently running operating system
    :return: system
    >>> get_system()
    'win'
    """
    import platform
    plat = platform.system()
    if plat.upper() == "WINDOWS":
        return "win"
    if plat.upper() == "DARWIN":
        return "mac"
    if plat.upper() == "LINUX":
        return "lin"


def post_update_things():
    """
    Performs some file operations after the update.
    :return:
    """
    import stat
    dsrdir = os.environ["DSR_DIR"]
    plat = get_system()
    upath = os.path.join(dsrdir, "dsr")
    try:
        if plat == "win":
            pass
        elif plat == "mac":
            shutil.copy2(os.path.abspath(os.path.join(dsrdir, "setup//dsr-mac")), upath)
            st = os.stat(upath)
            os.chmod(upath, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        else:
            shutil.copy2(os.path.abspath(os.path.join(dsrdir, "setup//dsr-linux")), upath)
            st = os.stat(upath)
            os.chmod(upath, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    except IOError:  # Unable to write in this case
        print('*** Unable to perform update. Please run me with super-user rights, e.g.: "sudo /opt/DSR/dsr -u" ***')


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
    try:
        dsrdir = os.path.dirname(os.path.realpath(__file__))
    except KeyError:
        print("*** Could not determine the location of DSR. Can not update. ***" )
        sys.exit()
    response = urllib.urlopen('{}/DSR-{}.tar.gz'.format(urlprefix, version))
    with tempfile.NamedTemporaryFile(delete=False) as tmpfile:
        tmpfile.write(response.read())
    tmpdir = tempfile.mkdtemp()  # a temporary directory
    try:
        with tarfile.open(tmpfile.name) as tarobj:
            tarobj.extractall(path=tmpdir)
    except tarfile.ReadError as e:
        print('*** Cound not get update from server. If this problem persists, please update manually! ***')
        print('***', e, '***')
        return False
    os.remove(tmpfile.name)
    try:
        overwrite_dir(os.path.join(tmpdir, "DSR-{}".format(version)), dsrdir, move=False)
    except OSError:
        print('*** Unable to perform update. Please run me with super-user rights, e.g.: "sudo /opt/DSR/dsr -u" ***')
        sys.exit()
    shutil.rmtree(tmpdir, ignore_errors=True)  # cleanup the files
    post_update_things()
    return True


if __name__ == "__main__":
    import doctest

    failed, attempted = doctest.testmod()  # verbose=True)
    if failed == 0:
        print('passed all {} tests!'.format(attempted))
    else:
        print('{} of {} tests failed'.format(failed, attempted))
