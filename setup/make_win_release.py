#!/usr/bin/python3
# -*- coding: utf-8 -*-

r"""
This script has to be run from the main dir e.g. D:\GitHub\DSR
"""
import os
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
    """
    Find the Inno Setup compiler (ISCC.exe), preferring Inno Setup 7 over 6.
    Exits with an error if Inno Setup is not installed.
    """
    uninstall_key = r'SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall'
    # Inno Setup 7 has a 64-bit and a 32-bit edition, both can be installed per machine or per user:
    app_ids = ['Inno Setup 7', 'Inno Setup 7 (32-bit)', 'Inno Setup 6']
    registry_views = [
        (winreg.HKEY_LOCAL_MACHINE, winreg.KEY_WOW64_64KEY),
        (winreg.HKEY_LOCAL_MACHINE, winreg.KEY_WOW64_32KEY),
        (winreg.HKEY_CURRENT_USER, 0),
    ]
    for app_id in app_ids:
        for hkey, view in registry_views:
            try:
                with winreg.OpenKey(hkey, rf'{uninstall_key}\{app_id}_is1', 0, winreg.KEY_READ | view) as key:
                    install_location, _ = winreg.QueryValueEx(key, 'InstallLocation')
            except OSError:
                continue
            compiler = pathlib.Path(install_location) / 'ISCC.exe'
            if compiler.exists():
                return str(compiler)
    # Fallback to common installation paths
    for version in ('7', '6'):
        for base_dir in (os.environ.get('ProgramFiles', r'C:\Program Files'),
                         os.environ.get('ProgramFiles(x86)', r'C:\Program Files (x86)'),
                         r'D:\Programme'):
            compiler = pathlib.Path(base_dir) / f'Inno Setup {version}' / 'ISCC.exe'
            if compiler.exists():
                return str(compiler)
    sys.exit('*** Inno Setup compiler (ISCC.exe) not found. Please install Inno Setup 7. ***')


def make_distribs():
    innosetup_compiler = get_innosetup_path()
    print(f'Using Inno Setup compiler: {innosetup_compiler}')
    # Run DSR setup compiler
    result = subprocess.run([innosetup_compiler, f'/dMyAppVersion={VERSION}', r'setup\dsr-install.iss', ])
    if result.returncode != 0:
        sys.exit(f'*** Inno Setup compiler failed with exit code {result.returncode}. ***')


# Make a zip file for web interface distribution:
# make_zip(files)

# Make binary distributions:
make_distribs()
