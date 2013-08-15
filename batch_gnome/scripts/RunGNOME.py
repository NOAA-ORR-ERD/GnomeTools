#!/usr/bin/env python

"""

RunGNOME

A script to run GNOME

NOTE: only works on Windows and with the GUI gnome on Mac at the moment.

"""

import sys, os, time, math
from TAP_Setup import setup

if not os.path.exists(setup.GNOME):
    raise Exception("The GNOME executable must be in the location specified in TAP_Setup (%s)"%setup.GNOME)

startTime = time.time()

##Start up GNOMES (in background)
for machine in range(1, setup.NumTrajMachines+1):
    CommandFilePath = os.path.join(setup.RootDir,'Machine' + str(machine), 'command.txt')
    CommandFilePathRel = os.path.relpath(CommandFilePath)
    if sys.platform == 'darwin':
        cmd = setup.GNOME + '/Contents/MacOS/gnome ' + CommandFilePathRel + " &"
        print cmd
        os.system(cmd)
    else: # assume Windows
        os.system('start ' + setup.GNOME + ' ' + CommandFilePathRel)
    time.sleep(5)

## check if they are all done running before exiting script
Trajectories = file(os.path.join(setup.RootDir,'ExpectedTrajectories.txt')).readlines()

#
##check for finished trajectories in order they are being done (i.e. by different machines)
#order = []
#RunsPerMachine = int(math.ceil(float(setup.NumStarts) / setup.NumTrajMachines))
#for m in range(0, RunsPerMachine):
#    order.extend(range(m, setup.NumStarts, RunsPerMachine))
#Trajectories = [ Trajectories[i] for i in order ]
#
#if setup.NumTrajMachines < setup.NumStarts:
#    ii = setup.NumTrajMachines #start from second file so we can get trajectory file size
#    FirstTime = 1
#    print 'Trajectory progress...'
#    while ii < len(Trajectories):
#        TrajectoryFile = Trajectories[ii].strip()
#        if os.path.exists(TrajectoryFile): #file exists
#            try:
#                if os.path.getmtime(TrajectoryFile) > startTime: #file is new
#                    if FirstTime and ii < setup.NumTrajMachines * 2: #first time
#                        expected_fsize = os.path.getsize(Trajectories[0].strip())
#                        print 'Done:  ' + Trajectories[ii - setup.NumTrajMachines].strip()
#                        FirstTime = 0
#                    if abs(os.path.getsize(TrajectoryFile)-expected_fsize) < 1000: #file is big enough
#                        print 'Done:  ' + TrajectoryFile.strip()
#                        ii = ii + 1
#                        FirstTime = 1
#                    else:
#                        time.sleep(3)
#            except WindowsError:
#                 time.sleep(3)
#else:
#    print 'Cube building may start before trajectories finish'
#    #Not sure how to check when done in this case -- unknown expected filesize


print "GNOMEs running: waiting for runs to be done"

# first wait for them all to exist
#Clean up whitespace
Trajectories = [f.strip() for f in Trajectories]
NotDoneTrajectories = Trajectories[:]
while NotDoneTrajectories:
    for filename in NotDoneTrajectories[:]:
        if os.path.exists(filename):
            if os.path.getmtime(filename) > startTime: #file is new
                NotDoneTrajectories.remove(filename)

print "All trajectory files exist: checking size"

# then check size
# assume the first one is done, and the right size
size = os.path.getsize(Trajectories[0])
AllSameSize = False
NotRightSizeTrajectories = NotDoneTrajectories = Trajectories[:]
while NotRightSizeTrajectories:
    for filename in NotRightSizeTrajectories[:]:
        if os.path.getsize(filename) == size:
            NotRightSizeTrajectories.remove(filename)
            
print "All trajectories are the right size"
                
                
                
        




