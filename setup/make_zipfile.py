"""
Creates a zip file with the content of the StructureDB program.
"""
from tarfile import TarFile

import os

from dsr import VERSION
from misc import copy_file, remove_file, walkdir

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
    os.chdir('../')
    zipfilen = 'setup/Output/DSR-{}.tar.gz'.format(version)
    remove_file(zipfilen)
    with TarFile(zipfilen, mode='w') as myzip:
        for f in filelist:
            print("Adding {}".format(f))
            for file in walkdir(f, exclude=['.pyc']):
                try:
                    myzip.add(file)
                except FileNotFoundError as e:
                    print(e)
                    print("#####################")
                    return False
    #copy_file(zipfilen, 'StructureFinder/scripts/Output/')
    print("File written to {}".format(zipfilen))

if __name__ == "__main__":
    make_zip(files)