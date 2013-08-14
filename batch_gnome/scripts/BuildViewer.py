#!/usr/bin/env python

"""
A simple script that copies all the cubes and everything into the right places

This has not been well set up to work universally. It's only been tested on one setup

"""
import os, shutil
from TAP_Setup import setup
import BatchGnome

TAPViewerDir = os.path.join(setup.RootDir, setup.TAPViewerPath)

#Check if TAP Viewer Dir exists:
if not os.path.isdir(TAPViewerDir):
    print "making new TAP Viewer Directory"
    os.mkdir(TAPViewerDir)

# copy the exe
shutil.copy(os.path.join(setup.TAPViewerSource, "TAP.exe"), TAPViewerDir)

# Check for TAPDATA
TAPDATADir = os.path.join(TAPViewerDir,"TAPDATA")
if not os.path.isdir(TAPDATADir):
    print "Making TAPDATA Directory"
    os.mkdir(TAPDATADir)

# copy the TAPCONFIG file
shutil.copy(os.path.join(setup.TAPViewerSource, "TAPCONFIG.txt"), TAPDATADir)

# copy the site.txt file
shutil.copy(os.path.join(setup.RootDir,"site.txt"), TAPDATADir)
shutil.copy(os.path.join(setup.TAPViewerSource, setup.MapFileName), TAPDATADir)

# copy the start times file (not required, but it's good to have it there
print setup.StartTimeFiles
for (filename, _) in setup.StartTimeFiles:
    shutil.copy(filename, TAPDATADir)

FullCubesPath = os.path.join(setup.RootDir, setup.CubesPath) 

for (season, junk) in setup.Seasons:
    SeasonPath = os.path.join(TAPDATADir,season)
    if not os.path.isdir(SeasonPath):
        print "Creating:", SeasonPath
        os.mkdir(SeasonPath)
    SeasonCubesPath = os.path.join(FullCubesPath,season)
    for name in os.listdir(SeasonCubesPath):
        print "Moving:", name
        shutil.move(os.path.join(SeasonCubesPath,name),
                     os.path.join(SeasonPath,name) )

