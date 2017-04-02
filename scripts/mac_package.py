import os
import subprocess
import time
import sys
import dsr
import misc
import selfupdate
try:  # Python2:
    import urllib2
    http_error = urllib2.HTTPError
except(ImportError, AttributeError):  # Python3:
    import urllib
    http_error = IOError

volname = "DSR-"+dsr.VERSION
dmgname = os.path.abspath("../setup/Output/{}.dmg".format(volname))
skeldmg = os.path.abspath("../setup/Output/DSR-skel-rw.dmg")
finaltmpdmg = os.path.abspath("../setup/Output/DSR-tmp-rw.dmg")

# Make copy of skeleton dmg file:
misc.copy_file(skeldmg, finaltmpdmg)

mountcommmand = ["hdiutil", "attach", "-readwrite", "-noverify", "-noautoopen", finaltmpdmg]
unmountcommmand = ["hdiutil", "detach", "/Volumes/DSR-{}".format(dsr.VERSION)]
convert_to_compressed = ["hdiutil", "convert", finaltmpdmg, "-format", "UDZO",
                         "-imagekey", "zlib-level=9", "-o", dmgname]

# Mount the image, unmount first, just in case:
subprocess.Popen(unmountcommmand, stderr=subprocess.PIPE)

print('Mounting .dmg file for mac deployment')
print(finaltmpdmg+"\n")
subprocess.call(mountcommmand)

# Download .tar.gz package from web server and extract to mounted volume:
try:
    selfupdate.get_update_package(version=dsr.VERSION, destdir='/Volumes/DSR-install/DSR', post=False)
except http_error:
    subprocess.call(["hdiutil", "detach", "/volumes/DSR-install"])
    print("Version {} is not present on xs3-data!".format(dsr.VERSION))
    # Clean temporary image:
    misc.remove_file(finaltmpdmg)
    sys.exit()

#########################################################
#  Do modification stuff here:
misc.copy_file("/Volumes/DSR-install/DSR/setup/dsr-mac", "/Volumes/DSR-install/DSR/dsr")
subprocess.call(["chmod", "755", "/Volumes/DSR-install/DSR/dsr"])

# Rename the volume:
renamecommand = ["diskutil", "rename", "DSR-install", "DSR-{}".format(dsr.VERSION)]
subprocess.call(renamecommand)
#########################################################

subprocess.call(unmountcommmand)
time.sleep(1)

# convert to compressed image:
misc.remove_file(dmgname)
subprocess.call(convert_to_compressed)

# Clean temporary image:
misc.remove_file(finaltmpdmg)

