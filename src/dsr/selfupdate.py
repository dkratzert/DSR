# -*- encoding: utf-8 -*-
# möp
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <dkratzert@gmx.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
from __future__ import print_function

import hashlib
import os
import shutil
import sys
import tarfile
import tempfile

from dsr import VERSION

urlprefix = "https://dkratzert.de/files/dsr"


from urllib.request import FancyURLopener


class MyOpener(FancyURLopener):
    """
    Sets the user agent of the urllib http request.
    """
    version = 'DSR cmdline {}'.format(VERSION)


myurlopen = MyOpener()


def get_current_dsr_version(silent=False):
    """
    determines the current version of DSR on the web server

    >>> get_current_dsr_version()
    '202'

    Returns
    -------
    version number
    :type: str
    """
    import socket
    socket.setdefaulttimeout(3)
    FancyURLopener.version = "DSR-updater {}".format(VERSION)
    try:
        response = myurlopen.open('{}/version.txt'.format(urlprefix))
    except IOError:
        if not silent:
            print("*** Unable to connect to update server. No Update possible. ***")
        return 0
    version = response.readline().decode('ascii').strip()
    return version


def is_update_needed(silent=False):
    # type: (bool) -> bool
    """
    Decides if an update of DSR is needed.
    :return: bool
    >>> is_update_needed()
    False
    """
    try:
        version = int(get_current_dsr_version(silent))
    except ValueError:
        print('*** Could not get version number from update server. ***')
        return False
    if int(VERSION) < int(version):
        return True
    else:
        return False


def update_dsr(force=False, version=None):
    """
    Updates the running DSR to the current version on the web server.
    #>>> update_dsr(force=True)
    #True
    """
    if not version:
        version = get_current_dsr_version()
    if not version:
        print('*** Could not get current version from server. ***')
    if force:
        status = get_update_package(version)
        if status:
            print('*** Finished updating to version {} ***'.format(version))
            return True
        else:
            return False
    if int(VERSION) < int(version):
        print('*** Current available version of DSR is {} you have installed version {}. '
              'Performing upate ... ***'.format(version, VERSION))
        status = get_update_package(version)
        if status:
            print('*** Finished updating to version {} ***'.format(version))
            return True
        else:
            print('*** Could not update DSR. ***')
            return False
    if (int(VERSION) >= int(version)) and int(version) > 0:
        print('*** DSR is already up to date (version {}) ***'.format(VERSION))
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


def post_update_things(dsrdir):
    """
    Performs some file operations after the update.
    :return: None
    """
    import stat
    plat = get_system()
    upath = os.path.join(dsrdir, "../../dsr")
    try:
        if plat == "win":
            pass
        elif plat == "mac":
            shutil.copy2(os.path.abspath(os.path.join(dsrdir, "../../setup/dsr-mac")), upath)
            st = os.stat(upath)
            os.chmod(upath, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        elif plat == "lin":
            shutil.copy2(os.path.abspath(os.path.join(dsrdir, "../../setup/dsr-linux")), upath)
            st = os.stat(upath)
            os.chmod(upath, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    except IOError:  # Unable to write in this case
        print('*** Unable to perform update. Please run me with super-user rights, e.g.: "sudo /opt/DSR/dsr -u" ***')


def move_dir(root_src_dir, root_dst_dir, move=True):
    """
    Moves the content of scrdir over destdir and overwrites all files.
    
    :param move: move directory instead of copying it
    :param root_src_dir: source directory
    :param root_dst_dir: target directory
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


def get_update_package(version, destdir=None, post=True):
    """
    Downloads the current DSR distribution from the web server and
    updates the files.

    :param post: Defines if post_update_things() should be executed
    :param destdir: Optional destdir instead of DSR_DIR
    :type version: int or string
    :return True/False

    >>> get_update_package('202')
    True
    """
    try:
        dsrdir = os.environ["DSR_DIR"]
    except KeyError:
        print("*** Could not determine the location of DSR. Can not update. ***")
        sys.exit()
    # DSR file:
    response = myurlopen.open('{}/DSR-{}.tar.gz'.format(urlprefix, version))
    try:
        with tempfile.NamedTemporaryFile(delete=False) as tmpfile:
            tmpfile.write(response.read())
    except Exception as e:
        print('*** Update timed out. Try again later. ***')
        print(e)
        return False
    downloaded_sha, tgz_sha = check_checksum(tmpfile, version)
    if downloaded_sha != tgz_sha:
        print('*** Checksum mismatch. Unable to update. If this problem persists, please update manually! ***')
        return False
    else:
        print("*** Checksums matched ***")
    # a temporary directory:
    tmpdir = tempfile.mkdtemp()
    try:
        with tarfile.open(tmpfile.name) as tarobj:
            tarobj.extractall(path=tmpdir)
    except tarfile.ReadError as e:
        print('*** Could not get update from server. If this problem persists, please update manually! ***')
        print('***', e, '***')
        return False
    os.remove(tmpfile.name)
    try:
        if not destdir:
            move_dir(os.path.join(tmpdir, "DSR-{}".format(version)), dsrdir, move=False)
        else:
            move_dir(os.path.join(tmpdir, "DSR-{}".format(version)), destdir, move=False)
    except OSError:
        print('*** Unable to perform update. Please run me with super-user rights, e.g.: "sudo /opt/DSR/dsr -u" ***')
        sys.exit()
    shutil.rmtree(tmpdir, ignore_errors=True)  # cleanup the files
    if post:
        post_update_things(dsrdir)
    return True


def check_checksum(tmpfile, version):
    """
    Creates a SHA512 checksum for the downloaded update package.
    :param tmpfile: dowloadad update package
    :param version: version numberof update
    :return: the two checksums
    """
    # download SHA file:
    response2 = myurlopen.open('{}/DSR-{}-sha512.sha'.format(urlprefix, version))
    downloaded_sha = response2.read()
    # Checksum for program package:
    tgz_sha = sha512_checksum(tmpfile.name)
    return downloaded_sha, tgz_sha


def sha512_checksum(filename, block_size=65536):
    """
    Calculates a SHA512 checksum from a file.

    :param filename:
    :param block_size:
    :return: str

    >>> sha512_checksum("../DSR-207.tar.gz")
    'e8d14033578e0ecce0d6c123a947060f9883fa735d1d3226b4f03f08a7eacecd'
    """
    sha512 = hashlib.sha512()
    with open(filename, 'rb') as f:
        for block in iter(lambda: f.read(block_size), b''):
            sha512.update(block)
    return sha512.hexdigest()


if __name__ == '__main__':
    import doctest

    failed, attempted = doctest.testmod()  # verbose=True)
    if failed == 0:
        print('passed all {} tests!'.format(attempted))
    else:
        print('{} of {} tests failed'.format(failed, attempted))
