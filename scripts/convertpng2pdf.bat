@echo off
mkdir output
C:\Python27\python "bildconvert.py" %1




rem #cd %1
rem #C:\Python27\python "..\convertpng2pdf.py"
rem cd ..

rem old crap:
rem for %%i in (*.png) do f:\ImageMagick-6.8.9-4\convert.exe %%i -label '%t' -geometry 600x600-30-40 test.png
rem f:\ImageMagick-6.8.9-4\montage.exe  -geometry 600x600 water.png -label '%t'   test.png
rem f:\ImageMagick-6.8.9-4\montage.exe *.png -label '%t' -tile 4x6 alle_bilder.pdf