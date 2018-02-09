"""
Creates a zip file with the content of the StructureDB program.
"""
import tarfile
import tempfile
import os
import shutil

from dsr import VERSION
from misc import copy_file, remove_file, walkdir
from selfupdate import sha256_checksum

version = VERSION

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
    fulldir = os.path.abspath(os.path.join(tmpdir, maindir))
    os.makedirs(fulldir)
    # Tar output file
    zipfilename = os.path.abspath('setup/Output/DSR-{}.tar.gz'.format(version))
    remove_file(zipfilename)
    # Go through DSR path and add files from list:
    for f in filelist:
        # Also add recoursive dirs:
        for filen in walkdir(f, exclude=['.pyc']):
            print(filen)
            # need path without filename to create target directories:
            path, _ = os.path.split(filen)
            target_dir = os.path.join(fulldir, path)
            if not os.path.exists(target_dir):
                os.makedirs(target_dir)
            copy_file(filen, target_dir)
    with tarfile.open(zipfilename, mode='w:gz') as archive:
        archive.add(fulldir, arcname=maindir, recursive=True)
    #copy_file(zipfilen, 'StructureFinder/scripts/Output/')
    print(fulldir)
    print("File written to {}".format(zipfilename))
    shutil.rmtree(tmpdir)
    make_shasum(zipfilename)

def make_shasum(filename):
    sha = sha256_checksum(filename)
    shafile = os.path.abspath('setup/Output/DSR-sha256-v{}.sha'.format(VERSION))
    with open(shafile, 'w') as f:
        f.write(sha)
    print("SHA256: {}".format(sha))

if __name__ == "__main__":
    make_zip(files)