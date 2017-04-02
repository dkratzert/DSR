import os
import subprocess
import time
import dsr
import misc
import selfupdate

volname = "DSR-"+dsr.VERSION
dmgname = os.path.abspath("../setup/Output/{}.dmg".format(volname))
skeldmg = os.path.abspath("../setup/Output/DSR-skel-rw.dmg")
finaltmpdmg = os.path.abspath("../setup/Output/DSR-tmp-rw.dmg")

print('Mounting .dmg file for mac deployment')
print(finaltmpdmg+"\n")

# Make copy of skeleton dmg file:
misc.copy_file(skeldmg, finaltmpdmg)

mountcommmand = ["hdiutil", "attach", "-readwrite", "-noverify", "-noautoopen", finaltmpdmg]
renamecommand = ["diskutil", "rename", "DSR-install", "DSR-{}".format(dsr.VERSION)]
unmountcommmand = ["hdiutil", "detach", "/volumes/DSR-{}".format(dsr.VERSION)]
convertfinal = ["hdiutil", "convert", finaltmpdmg, "-format", "UDZO", "-imagekey", "zlib-level=9", "-o", dmgname]

# Mount the image, unmount first, just in case:
subprocess.call(unmountcommmand)
p = subprocess.Popen(mountcommmand, stdout=subprocess.PIPE)
line = p.stdout.readlines()

# Download .tar.gz package from web server and extract to mounted volume:
selfupdate.get_update_package(version=dsr.VERSION, destdir='/volumes/DSR-install/DSR', post=False)

# Do modification stuff here:
misc.copy_file("/volumes/DSR-install/DSR/setup/dsr-mac", "/volumes/DSR-install/DSR/dsr")
subprocess.call(["chmod", "755", "/volumes/DSR-install/DSR/dsr"])

# Rename the volume:
subprocess.call(renamecommand)

subprocess.call(unmountcommmand)
time.sleep(1)

# convert to compressed image:
subprocess.call(convertfinal)

# Clean temporary image:
misc.remove_file(finaltmpdmg)

