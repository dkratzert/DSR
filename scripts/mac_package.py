import os
import subprocess
import tempfile
import misc
import selfupdate

try:  # Python2:
    # noinspection PyCompatibility
    import urllib2
    http_error = urllib2.HTTPError
except(ImportError, AttributeError):  # Python3:
    import urllib
    http_error = IOError

version = "207"  # define the .dmg release version

volname = "DSR-"+version
inputfile = "../setup/Output/DSR-{}.tar.gz".format(version, version)
dmgname = os.path.abspath("../setup/Output/{}.dmg".format(volname))
skeldmg = os.path.abspath("../setup/Output/DSR-skel-rw.dmg")
finaltmpdmg = os.path.abspath("../setup/Output/DSR-tmp-rw.dmg")

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

"""
# use instead of direct file opening below:
# Download .tar.gz package from web server and extract to mounted volume:
try:
    print("getting DSR from web server ....")
    selfupdate.get_update_package(version=version, destdir='/Volumes/DSR-install/DSR', post=False)
except http_error:
    subprocess.call(["hdiutil", "detach", "/volumes/DSR-install"])
    print("Version {} is not present on xs3-data!".format(version))
    # Clean temporary image:
    misc.remove_file(finaltmpdmg)
    sys.exit()
"""
# use instead of webserver stuff above
tmpdir = tempfile.mkdtemp()  # a temporary directory
misc.extract_tarfile(inputfile, tmpdir)
selfupdate.move_dir(os.path.join(tmpdir, "DSR-{}".format(version)), '/Volumes/DSR-install/DSR/')
# end

#########################################################
#  Do modification stuff here:
misc.copy_file("/Volumes/DSR-install/DSR/setup/dsr-mac", "/Volumes/DSR-install/DSR/dsr")
subprocess.call(["chmod", "755", "/Volumes/DSR-install/DSR/dsr"])

# Rename the volume:
renamecommand = ["diskutil", "rename", "DSR-install", "DSR-{}".format(version)]
subprocess.call(renamecommand)
#########################################################

subprocess.call(unmountcommmand)
#time.sleep(1)

# convert to compressed image:
misc.remove_file(dmgname)
subprocess.call(convert_to_compressed)

# Clean temporary image:
misc.remove_file(finaltmpdmg)

