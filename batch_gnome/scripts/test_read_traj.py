"""
module to test trajectory reading
"""

import numpy as np

print "running test"
import TAP_mod
reload(TAP_mod)

import TAP_Setup as setup

print "run length: %s hours, %s days"%(setup.TrajectoryRunLength, setup.TrajectoryRunLength/24.)

traj_file = r"..\..\Trajectories\Spring\001\time001.bin"
#traj_file = r"..\..\Trajectories\test2.bin"

(Trajectory,(NumTimesteps,NumLEs),HeaderData,flags) = TAP_mod.ReadTrajectory(traj_file)
    
print  "NumTimesteps", NumTimesteps
print "NumLEs", NumLEs
print "HeaderData", HeaderData

print "traj shape", Trajectory.shape
print "path of LE 0:"
#print Trajectory[:,0,:5]
print "flags", flags[:,1]

## raw read:
d = file(traj_file, 'rb').read()
print repr(d[400:467])
ind = d.find("DATA]\r\n") + 7 


rec = d[ind:ind+17]
print "first record"
print np.fromstring(rec[:4], dtype=np.int32)
print np.fromstring(rec[4:8], dtype=np.int32)
print np.fromstring(rec[8:12], dtype=np.float32)
print np.fromstring(rec[12:13], dtype=np.uint8)
rec = d[ind+17:ind+17+17]
print "second record"
print np.fromstring(rec[:4], dtype=np.int32)
print np.fromstring(rec[4:8], dtype=np.int32)
print np.fromstring(rec[8:12], dtype=np.float32)
print np.fromstring(rec[12:13], dtype=np.uint8)

print "long:", Trajectory[:,:,0]
print "lat:", Trajectory[:,:,1]
