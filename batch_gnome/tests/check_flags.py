#!/usr/bin/env python

"""
check flags in a TAP binary file
"""

import sys

from batch_gnome import tap_mod

infilename = sys.argv[1]

(Trajectory,(NumTimesteps,NumLEs),HeaderData,flags) = tap_mod.ReadTrajectory(infilename)

print "flags are:"
print flags

