#!/usr/bin/env python

filename = "Data/2005110317_CHL_a.gbinwinds_42036_d300_h36_d300_h36.cts_0.04.bout"

import TAP_mod
import numpy as N

(Trajectory,(NumTimesteps,NumLEs),HeaderData,flags) = TAP_mod.ReadTrajectory(filename, CheckFlag = 1)
(Trajectory2,(NumTimesteps,NumLEs),HeaderData,flags) = TAP_mod.ReadTrajectory_Old(filename, CheckFlag = 1)

print HeaderData


if N.alltrue(Trajectory == Trajectory2):
    print "They are the same!"
else:
    print "OOPS! -- something is wrong"

