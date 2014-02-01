#/usr/bin/env python
#-*- encoding: utf-8 -*-


from py2exe.build_exe import py2exe 
from distutils.core import setup

#import sys,os
"""
origIsSystemDLL = py2exe.build_exe.isSystemDLL()
def isSystemDLL(pathname):
        if os.path.basename(pathname).lower() in ("USER32.dll", "SHELL32.dll" ,"ADVAPI32.dll", "WS2_32.dll", "GDI32.dll", "KERNEL32.dll"):
                return 0
        return origIsSystemDLL(pathname)
py2exe.build_exe.isSystemDLL = isSystemDLL
"""

setup( 
      console=["dsr.py", "atoms.py"] 
      )
