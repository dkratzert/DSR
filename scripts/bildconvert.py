#-*- encoding: utf-8 -*-

'''
Anschließend:
montage.exe *.png  -geometry 600x600+1+1 -tile 3x6 test.pdf

montage.exe *.png -trim +repage -geometry 600x600+10+10 -tile 4x6 test.pdf

montage.exe *.png -geometry 600x600-30-40 -tile 4x6 test.pdf
'''

import os, sys
from subprocess import call
import fnmatch

try:
    outdir = sys.argv[1]
except(IndexError):
    print('Please give the output dir as argument.')
    sys.exit()

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


prog = "f:/programme/Java/bin/java"

files = os.listdir('.')

list_of_fragments = []


for f in files:
    if fnmatch.fnmatch(f, '*.res'):
        list_of_fragments.append(f)

#list_of_fragments = ['toluene']

def set_options(fragment):
    options = ' -Djmol.logger.info=false -Djmol.logger.warn=false -Xmx512m -jar \
        F:\Programme\jmol\Jmol.jar {}.res -ions bildererstellung.txt -g 600x600 \
        -w PNG:{}\{}.png -x'.format(fragment, outdir, fragment)
    #options = ' -Xmx512m  -jar "F:\Programme\jmol\jmol.jar"'
    #options = ""
    return options

for num, fragment in enumerate(list_of_fragments):
    fileName, fileExtension = os.path.splitext(fragment)
    options = set_options(fileName)
    num = str(num+1)
    print('\n'+fileName+'  '+num+'\n')
    os.system(prog+options)


