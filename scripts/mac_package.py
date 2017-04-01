import os
import subprocess
import sys

import time

import dsr
import misc

try:
    src = os.path.abspath(sys.argv[1])
except IndexError:
    print("You have to provide the path to the deploy directory as first argument")
    sys.exit()

volname = "DSR-"+dsr.VERSION
dmgname = volname+".dmg"
dmgtmp = volname+"-tmp.dmg"
dsrbat = os.path.join(src, "DSR", "dsr")
print(dsrbat)

# use tar.gz file instead of directory
# use DSR/setup/output as target directory

misc.copy_file(os.path.join(src, "DSR", "setup", "dsr-mac"), dsrbat)
subprocess.call(["chmod", "755", dsrbat])

misc.remove_file(dmgtmp)

command_line = ["hdiutil", "create", dmgtmp, "-srcfolder", src, "-fs", "HFS+",
                "-format", "UDRW", "-size", "20M", "-volname", volname]
# Maybejust use the DSR-skel instead of create:
print('Creating .dmg file for mac deployment')
print("\n"+' '.join(command_line)+"\n")
subprocess.call(command_line)
# mount the image:
mountcommmand = ["hdiutil", "attach", "-readwrite", "-noverify", "-noautoopen", dmgtmp]
unmountcommmand = ["hdiutil", "detach", "/volumes/{}".format(volname)]

subprocess.call(mountcommmand)

# do modification stuff here
time.sleep(1)
subprocess.call(unmountcommmand)
time.sleep(1)
# convert to compressed image:
convertfinal = ["hdiutil", "convert", dmgtmp, "-format", "UDZO", "-imagekey", "zlib-level=9", "-o", dmgname]
print("\n"+' '.join(convertfinal))
subprocess.call(convertfinal)
misc.remove_file(dmgtmp)

