mkdir %1

C:\Python27\python "bildconvert.py" %1

f:\ImageMagick-6.8.9-4\montage.exe %1\*.png -geometry 600x600-30-40 -tile 4x6 alle_bilder.pdf
