#!/usr/bin/env python

"""
test script for volume cube 
"""
import os, sys
import numpy as np


## fixme: this really needs a better structure
sys.path.append( os.path.split(os.path.split(os.path.abspath(__file__))[0])[0])

from lib import TAP_mod, tap_comp_volume
reload(tap_comp_volume)
reload(TAP_mod)


# some test data
infilename = "micro_test.bin"
                                                
(Trajectory,(NumTimesteps,NumLEs),HeaderData,Flags) = TAP_mod.ReadTrajectory(infilename)

print Flags
print Flags.dtype

max = Trajectory.max(axis=0).max(axis=0)
min = Trajectory.min(axis=0).min(axis=0)

# set up a grid:
grid = TAP_mod.Grid(min[0], max[0], min[1], max[1], 10, 10)
print grid

# loop through the time steps:
for t_ind in xrange(NumTimesteps):
#for t_ind in [10]:
    print "timestep :", t_ind
    Positions = Trajectory[t_ind]
    flags = Flags[t_ind]
    masses = np.ones_like(flags, dtype=np.float32)
    Volumes = tap_comp_volume.comp_volume(Positions, masses, flags, grid)
    print Volumes

