#!/usr/bin/python3
# -*- coding: utf-8 -*-

r"""
This script has to be run from the main dir e.g. D:\GitHub\StructureFinder
"""
import pathlib
import subprocess
import sys
import winreg

pth = pathlib.Path(__file__).parent.parent
print(pth)
sys.path.insert(0, str(pth / 'src/dsr_shelx'))
sys.path.insert(0, str(pth))

from version import VERSION

print("Updating version numbers to version {} ...".format(VERSION))

#print("Linux deb... {}".format(VERSION))
#process_debian_and_spec(debianpath)

#print("windows iss... {}".format(VERSION))
#process_iss(isspath)

print("Version numbers updated.")


def get_innosetup_path() -> str:
    """Find the Inno Setup compiler by querying the Windows registry."""
    registry_keys = [
        (winreg.HKEY_LOCAL_MACHINE, r'SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\Inno Setup 6_is1'),
        (winreg.HKEY_LOCAL_MACHINE, r'SOFTWARE\WOW6432Node\Microsoft\Windows\CurrentVersion\Uninstall\Inno Setup 6_is1'),
        (winreg.HKEY_CURRENT_USER, r'Software\Microsoft\Windows\CurrentVersion\Uninstall\Inno Setup 6_is1'),
    ]
    for hkey, subkey in registry_keys:
        try:
            with winreg.OpenKey(hkey, subkey) as key:
                install_location, _ = winreg.QueryValueEx(key, 'InstallLocation')
                compiler = pathlib.Path(install_location) / 'ISCC.exe'
                if compiler.exists():
                    return str(compiler)
        except OSError:
            continue
    # Fallback to common installation paths
    fallback_paths = [
        pathlib.Path(r'C:\Program Files (x86)\Inno Setup 6\ISCC.exe'),
        pathlib.Path(r'C:\Program Files\Inno Setup 6\ISCC.exe'),
        pathlib.Path(r'D:\Programme\Inno Setup 6\ISCC.exe'),
    ]
    for path in fallback_paths:
        if path.exists():
            return str(path)
    raise FileNotFoundError('Inno Setup 6 compiler (ISCC.exe) not found. Please install Inno Setup 6.')


def make_distribs():
    innosetup_compiler = get_innosetup_path()
    # Run DSR setup compiler
    subprocess.call([innosetup_compiler, f'/dMyAppVersion={VERSION}', r'setup\dsr-install.iss', ])


# Make a zip file for web interface distribution:
# make_zip(files)

# Make binary distributions:
make_distribs()
