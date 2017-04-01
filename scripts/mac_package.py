import os
import subprocess
import stat
import sys

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

misc.copy_file(os.path.join(src, "DSR", "setup", "dsr-mac"), dsrbat)
subprocess.call(["chmod", "755", dsrbat])

misc.remove_file(dmgtmp)

command_line = ["hdiutil", "create", dmgtmp, "-srcfolder", src, "-fs", "HFS+",
                "-format", "UDRW", "-volname", volname]

print('Creating .dmg file for mac deployment')
print("\n"+' '.join(command_line)+"\n")
subprocess.call(command_line)

mountcommmand = ["hdiutil", "attach", "-readwrite", "-noverify", "-noautoopen", dmgtmp]
subprocess.call(mountcommmand)

"""

After lots of research, I've come up with this answer, and I'm hereby putting it here as an answer for my own question, for reference:

Make sure that "Enable access for assistive devices" is checked in System Preferences>>Universal Access. It is required for the AppleScript to work. You may have to reboot after this change (it doesn't work otherwise on Mac OS X Server 10.4).
Create a R/W DMG. It must be larger than the result will be. In this example, the bash variable "size" contains the size in Kb and the contents of the folder in the "source" bash variable will be copied into the DMG:
hdiutil create -srcfolder "${source}" -volname "${title}" -fs HFS+ \
      -fsargs "-c c=64,a=16,e=16" -format UDRW -size ${size}k pack.temp.dmg
Mount the disk image, and store the device name (you might want to use sleep for a few seconds after this operation):
device=$(hdiutil attach -readwrite -noverify -noautoopen "pack.temp.dmg" | \
         egrep '^/dev/' | sed 1q | awk '{print $1}')
Store the background picture (in PNG format) in a folder called ".background" in the DMG, and store its name in the "backgroundPictureName" variable.
Use AppleScript to set the visual styles (name of .app must be in bash variable "applicationName", use variables for the other properties as needed):
echo '
   tell application "Finder"
     tell disk "'${title}'"
           open
           set current view of container window to icon view
           set toolbar visible of container window to false
           set statusbar visible of container window to false
           set the bounds of container window to {400, 100, 885, 430}
           set theViewOptions to the icon view options of container window
           set arrangement of theViewOptions to not arranged
           set icon size of theViewOptions to 72
           set background picture of theViewOptions to file ".background:'${backgroundPictureName}'"
           make new alias file at container window to POSIX file "/Applications" with properties {name:"Applications"}
           set position of item "'${applicationName}'" of container window to {100, 100}
           set position of item "Applications" of container window to {375, 100}
           update without registering applications
           delay 5
           close
     end tell
   end tell
' | osascript


chmod -Rf go-w /Volumes/"${title}"
sync
sync
hdiutil detach ${device}
hdiutil convert "/pack.temp.dmg" -format UDZO -imagekey zlib-level=9 -o "${finalDMGName}"
rm -f /pack.temp.dmg
"""