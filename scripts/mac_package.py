import os
import subprocess
import tempfile

from src.dsr.dsr import VERSION
from scripts.make_zipfile import files, make_zip
from src.dsr import misc, selfupdate

try:  # Python2:
    # noinspection PyCompatibility
    import urllib2

    http_error = urllib2.HTTPError
except(ImportError, AttributeError):  # Python3:
    import urllib

    http_error = IOError

version = VERSION

# First, create a .tar.gz file with the DSR program:
make_zip(files)

volname = "DSR-" + version
inputfile = "setup/Output/DSR-{}.tar.gz".format(version)
dmgname = os.path.abspath("setup/Output/{}.dmg".format(volname))
skeldmg = os.path.abspath("setup/Output/DSR-skel-rw.dmg")
finaltmpdmg = os.path.abspath("setup/Output/DSR-tmp-rw.dmg")

# Make copy of skeleton dmg file:
misc.copy_file(skeldmg, finaltmpdmg)

mountcommmand = ["hdiutil", "attach", "-readwrite", "-noverify", "-noautoopen", finaltmpdmg]
unmountcommmand = ["hdiutil", "detach", "/Volumes/DSR-{}".format(version)]
convert_to_compressed = ["hdiutil", "convert", finaltmpdmg, "-format", "UDZO",
                         "-imagekey", "zlib-level=9", "-o", dmgname]

# Mount the image, unmount first, just in case:
subprocess.Popen(unmountcommmand, stderr=subprocess.PIPE)

print('Mounting .dmg file for mac deployment')
print('Mounting: ', finaltmpdmg, '\n')
subprocess.call(mountcommmand)

# use instead of webserver stuff above
tmpdir = tempfile.mkdtemp()  # a temporary directory
misc.extract_tarfile(inputfile, tmpdir)
selfupdate.move_dir(os.path.join(tmpdir, "DSR-{}".format(version)), '/Volumes/DSR-install/DSR/')
# end

#########################################################
#  Do modification stuff here:
misc.copy_file("/Volumes/DSR-install/DSR/setup/dsr-mac", "/Volumes/DSR-install/DSR/dsr")
misc.remove_file("/Volumes/DSR-install/DSR/dsr.bat")
subprocess.call(["chmod", "755", "/Volumes/DSR-install/DSR/dsr"])

# Rename the volume:
renamecommand = ["diskutil", "rename", "DSR-install", "DSR-{}".format(version)]
subprocess.call(renamecommand)
#########################################################

subprocess.call(unmountcommmand)
# time.sleep(1)

# convert to compressed image:
misc.remove_file(dmgname)
subprocess.call(convert_to_compressed)

# Clean temporary image:
misc.remove_file(finaltmpdmg)
