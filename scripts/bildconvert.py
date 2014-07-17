#-*- encoding: utf-8 -*-

'''
f:/programme/Java/bin/java -Djmol.logger.info=false -Djmol.logger.warn=false -Xmx512m -jar F:\Programme\jmol\Jmol.jar pfanion.res -ions bildererstellung.txt -g 600x600 -w PNG:pfanion.png -x
'''

import os, sys
from subprocess import call
import fnmatch



#try:
#    noconvert = sys.argv[1]
#except(IndexError):
#    pass
#

def write_options_file(fragment, options):
    ## write to file:
    optfile = fragment+'.txt'
    try:
        f = open(optfile, 'w')
        for line in options:
            f.write(line)
    except(IOError):
        print('could not write file {}'.format(optfile))
        sys.exit(-1)
    f.close()




files = os.listdir('.')
files.sort()
list_of_fragments = []


for f in files:
    if fnmatch.fnmatch(f, '*.res'):
        list_of_fragments.append(f)

#list_of_fragments = ['toluene']

def set_options(fragment):
    options = ' -Djmol.logger.info=false -Djmol.logger.warn=false -Xmx512m -jar \
        F:\Programme\jmol\Jmol.jar {}.res -ions bildererstellung.txt -g 600x600 \
        -w PNG:{}.png -x'.format(fragment, fragment)
    #options = ' -Xmx512m  -jar "F:\Programme\jmol\jmol.jar"'
    #options = ""
    return options



for num, fragment in enumerate(list_of_fragments):
    Name, fileExtension = os.path.splitext(fragment)
    java_prog = "f:/programme/Java/bin/java"
    java_options = set_options(Name)
    options_png = " {}.png -geometry 600x600-30-40-30-30 -font Arial -pointsize 40 label:{} \
                        -gravity South -composite {}.png".format(Name, Name, Name)
    prog_png = "f:\ImageMagick-6.8.9-4\convert.exe"
    num = str(num+1)
    print('\n'+Name+'  '+num+'\n')
    # remove next lines if only montage is required
    # create png images from res file:
    os.system(java_prog+java_options)
    # label image:
    os.system(prog_png+options_png)



# make tiles
optionsm = " *.png -geometry 600x600 -tile 3x4 output\\alle_bilder.png"
progm = "f:\ImageMagick-6.8.9-4\montage.exe"

os.system(progm+optionsm)
