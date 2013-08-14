#!/usr/bin/env python

"""
check flags in a TAP binary file
"""

import sys

sys.path.append("../lib") # to get the BatchGNOME lib stuff
import TAP_mod

infilename = sys.argv[1]

(Trajectory,(NumTimesteps,NumLEs),HeaderData,flags) = TAP_mod.ReadTrajectory(infilename)

print "flags are:"
print flags