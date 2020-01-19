#! /usr/bin/env python2.7

##############################################################################
# The program to superimpose atoms of two molecules by quaternion method
#
# David J. Heisterberg
# The Ohio Supercomputer Center
# 1224 Kinnear Rd.
# Columbus, OH  43212-1163
# (614)292-6036
# djh@ccl.net    djh@ohstpy.bitnet    ohstpy::djh
#
# Translated to C from fitest.f program and interfaced with Xmol program
# by Jan Labanowski,  jkl@ccl.net   jkl@ohstpy.bitnet   ohstpy::jkl
#
# Translated to python from quatfit.c
# by Thomas J. L. Mustard, mustardt@onid.orst.edu
#
# Copyright: Ohio Supercomputer Center, David J. Heisterberg, 1990.
# The program can be copied and distributed freely, provided that
# this copyright in not removed. You may acknowledge the use of the
# program in published material as:
# David J. Heisterberg, 1990, unpublished results.

from __future__ import print_function
import getopt
import math
import os
import shutil
import sys

import quatfit

### --- Arguments --- ###
program = "QuatfitPy-CLASSICAL.py"
reffile = ''
fitfile = ''
ofile = ''
pairsfile = ''
statfile = ''
debug = 0

# Read command line args
try:
    myopts, args = getopt.getopt(sys.argv[1:], "r:f:p:o:s:hd")
except getopt.GetoptError:
    print(program + " -r <ref> -f <fit> -p <pairs> -o <out> -s <stat> -h -d")
    sys.exit(2)
###############################
# o == option
# a == argument passed to the o
###############################
for o, a in myopts:
    if o == '-r':
        reffile = a
    elif o == '-f':
        fitfile = a
    elif o == '-o':
        ofile = a
    elif o == '-p':
        pairsfile = a
    elif o == '-s':
        statfile = a
    elif o == '-d':
        debug += 1
    elif o == '-h':
        print(program + " -r <ref> -f <fit> -p <pairs> -o <out> -s <stat> -h -d")
        sys.exit(0)
    else:
        print("Usage: %s -r <ref> -f <fit> -p <pairs> -o <out> -s <stat> -h -d" % sys.argv[0])
        sys.exit(0)

reffilelol = quatfit.parseXYZ(reffile)
fitfilelol = quatfit.parseXYZ(fitfile)

pairs, weights = quatfit.parsePairs(pairsfile)

fitfilelol, rms = quatfit.quatfitGetMolecule(reffilelol, fitfilelol, pairs, weights)

print("Weighted root mean square = " + str(rms))

quatfit.outputXYZ(ofile, fitfilelol, 1)
