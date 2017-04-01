import os
import subprocess

import sys

import dsr
import misc

try:
    src = os.path.abspath(sys.argv[1])
except IndexError:
    print("You have to provide the path to the deploy directory as first argument")
    sys.exit()

volname = "DSR-"+dsr.VERSION

misc.copy_file(os.path.join(src, "setup", "dsr-mac"), os.path.join(src, "dsr"))

os.fchmod(os.path.join(src, "dsr"), os.stat.S_IXUSR | os.stat.S_IXGRP | os.stat.S_IXOTH)

command_line = ["hdiutil", "create", "-srcfolder", src, "-format", "UDZO", "-volname", volname]

print('Creating .dmg file for mac deployment')
subprocess.Popen(command_line)
