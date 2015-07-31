#!/usr/bin/env python
"""

BuildAll


A simple script that runs all the TAP scripts
"""

import os, sys

try:
    RootDir = sys.argv[1]
except IndexError:
    raise Exception("You must pass in the RootDir on the command line")

if not os.path.exists(RootDir):
    raise Exception("RootDir: %s Doesn't exist"%RootDir)


os.system("python BuildStartTimes.py %s"%RootDir)

print
print "Running GNOME"
#NOTE: RunGNOME script needs some fixing so that it can better tell
#      That it's done -- it may move on too soon now.
os.system("python RunPyGnome.py %s"%RootDir)
#i = raw_input("Run GNOME now, then hit <Enter> to continue, <ctrl+C> to stop")

print 
print "Building the Cubes"
os.system("python BuildCubes.py %s"%RootDir)

print 
print "Building the Site.txt file"
os.system("python BuildSite.txt.py %s"%RootDir)

print 
print "Setting up the TAP viewer"
os.system("python BuildViewer.py %s"%RootDir)

print "Done"
