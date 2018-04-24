"""
Creates a zip file with the content of the StructureDB program.
"""
from __future__ import print_function
import string
import tarfile
import tempfile
import os
import shutil
import sys

from dsr import VERSION
from misc import copy_file, remove_file, walkdir
from selfupdate import sha512_checksum

files = [
    "afix.py",
    "README",
    "atomhandling.py",
    "atoms.py",
    "constants.py",
    "dbfile.py",
    "dsr.py",
    "dsrparse.py",
    "export.py",
    "misc.py",
    "elements.py",
    "options.py",
    "resfile.py",
    "restraints.py",
    "cf3fit.py",
    "selfupdate.py",
    "terminalsize.py",
    "resi.py",
    "refine.py",
    "pyperclip.py",
    "dsr_db.txt",
    "manuals/DSR-manual.pdf",
    "setup/dsr.sh",
    "setup/dsr-mac",
    "setup/dsr-linux",
    "dsr",
    "example/p21c.hkl",
    "example/p21c.res",
    "example/p21c_step0.res",
    "example/p21c_step1.res",
    "example/p21c_step2.ins",
    "example/p21c_step3.res",
    "example/p21c_final.res",
    "example/p21n_cf3.hkl",
    "example/p21n_cf3.res",
    "networkx",
    "mpmath"
]


def make_zip(filelist):
    """
    :type filelist: list
    """
    # Main path in tarfile:
    maindir = 'DSR-{}/'.format(VERSION)
    # tmpdir for file collection:
    tmpdir = tempfile.mkdtemp()
    # full directory in temp:
    full_tmp_dir = os.path.abspath(os.path.join(tmpdir, maindir))
    os.makedirs(full_tmp_dir)
    # Tar output file
    zipfilename = os.path.abspath('setup/Output/DSR-{}.tar.gz'.format(VERSION))
    remove_file(zipfilename)
    # Go through DSR path and add files from list:
    for f in filelist:
        # Also add recoursive dirs:
        for filen in walkdir(f, exclude=['.pyc']):
            print(filen)
            # need path without filename to create target directories:
            path, _ = os.path.split(filen)
            target_tmp_dir = os.path.join(full_tmp_dir, path)
            if not os.path.exists(target_tmp_dir):
                os.makedirs(target_tmp_dir)
            copy_file(filen, target_tmp_dir)
            dos2unix(os.path.join(target_tmp_dir, os.path.split(filen)[1]))  # dos2unix only in target tmp
    with tarfile.open(zipfilename, mode='w:gz') as archive:
        archive.add(full_tmp_dir, arcname=maindir, recursive=True)
    # copy_file(zipfilen, 'StructureFinder/scripts/Output/')
    print("\nFile written to {}".format(zipfilename))
    shutil.rmtree(tmpdir)
    make_shasum(zipfilename)


def make_shasum(filename):
    sha = sha512_checksum(filename)
    shafile = os.path.abspath('setup/Output/DSR-{}-sha512.sha'.format(VERSION))
    remove_file(shafile)
    with open(shafile, 'w') as f:
        f.write(sha)
    print("SHA512: {}".format(sha))


def dos2unix(filename):
    """
    >>> dos2unix('./profiling.bat')
    """
    if is_binary(filename):
        print('Binary file {} ignored.'.format(filename))
        return
    if sys.version_info[0] > 2:
        file_contents = open(filename, "r").read()
        f = open(filename, "w")
        f.write(file_contents)
        f.close()
    else:
        text = open(filename, 'rb').read().replace(b'\r\n', b'\n')
        open(filename, 'wb').write(text)


def istext(filename):
    """
    Got this from https://stackoverflow.com/questions/1446549/
                    how-to-identify-binary-and-text-files-using-python
    Not compatible to Python3
    """
    s=open(filename).read(1024)
    text_characters = "".join(map(chr, range(32, 127)) + list("\n\r\t\b"))
    _null_trans = string.maketrans("", "")
    if not s:
        # Empty files are considered text
        return True
    if "\0" in s:
        # Files with null bytes are likely binary
        return False
    # Get the non-text characters (maps a character to itself then
    # use the 'remove' option to get rid of the text characters.)
    t = s.translate(_null_trans, text_characters)
    # If more than 30% non-text characters, then
    # this is considered a binary file
    if float(len(t))/float(len(s)) > 0.30:
        return False
    return True


def is_binary(filename):
    """
    Return true if the given filename is binary.

    @raise EnvironmentError: if the file does not exist or cannot be accessed.
    @attention: found @ http://bytes.com/topic/python/answers/21222-determine-file-type-binary-text on 6/08/2010
    @author: Trent Mick <TrentM@ActiveState.com>
    @author: Jorge Orpinel <jorge@orpinel.com>"""
    fin = open(filename, 'rb')
    try:
        CHUNKSIZE = 1024
        while 1:
            chunk = fin.read(CHUNKSIZE)
            if b'\0' in chunk: # found null byte
                return True
            if len(chunk) < CHUNKSIZE:
                break # done
    # A-wooo! Mira, python no necesita el "except:". Achis... Que listo es.
    finally:
        fin.close()
    return False


if __name__ == "__main__":
    make_zip(files)
