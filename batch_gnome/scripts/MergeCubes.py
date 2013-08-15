#!/usr/bin/env python

import sys
import numpy as np


cube1name = sys.argv[1]
cube2name = sys.argv[2]
outfilename = sys.argv[3]

shape = (NumTimes, NumSites, NumSpills) = 12, 2046, -1

print "Using hard-coded shape: (%i, %i, %i)"%shape

print "merging: %s and %s, creating %s"%(cube1name, cube2name, outfilename)

outfile = file(outfilename, 'wb')

print "reading the cubes"
cube1 = np.fromfile(cube1name, dtype=np.uint8).reshape(shape)
cube2 = np.fromfile(cube2name, dtype=np.uint8).reshape(shape)

print "merging:"
cube3 = np.concatenate((cube1, cube2), axis=2)

print "writing:"
cube3.tofile(outfilename)

