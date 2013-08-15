#!/usr/bin/env python


import os, string, math
from datetime import datetime, timedelta

from TAP_Setup import setup

from batch_gnome import batch_gnome

# create command file dir if it doesn't exist
print "RootDir is", setup.RootDir

print "Setting up Master Command File"
cfile = batch_gnome.CommandFile()
cfile.SaveFilePath   = os.curdir
cfile.SaveFileName   = setup.SaveFileName
cfile.RunLength      = setup.TrajectoryRunLength
cfile.NumLEs         = setup.NumLEs
cfile.ModelTimeStep  = setup.ModelTimeStep

NumCubes = len(setup.CubeStartSites)
NumStarts = setup.NumStarts
NumSeasons = len(setup.StartTimeFiles)

RunsPerMachine = int( math.ceil(float((NumCubes * NumSeasons * NumStarts)) / setup.NumTrajMachines ) )
print "RunsPerMachine", RunsPerMachine
    
#CubesPerMachine = int( math.ceil(float(len(setup.CubeStartSites))  * NumSeasons / setup.NumTrajMachines) )

TrajectoryFile = file(os.path.join(setup.RootDir, 'ExpectedTrajectories.txt'),'w')
MachineNum = 0
RunNum = 0
for Season in setup.StartTimeFiles:
    SeasonName = Season[1]
    StartTimes = open(Season[0],'r').readlines()[:setup.NumStarts]
    CubeNum = 0
    for Site in setup.CubeStartSites:
        CubeNum += 1
        OutputFilePath  = os.path.abspath(os.path.join(setup.RootDir, setup.TrajectoriesPath, SeasonName, str(CubeNum).zfill(3) ))
        #print "OutputFilePath", OutputFilePath
        #print os.path.abspath(OutputFilePath)
        #raise Exception("stopping here")
        for i,time in zip(range(len(StartTimes)),StartTimes):
            if (RunNum % RunsPerMachine) == 0: # make a new command file
                MachineNum += 1
                if cfile.Runs:
                    cfile.write()
                    cfile.ClearRuns()
                MachinePath = os.path.join(setup.RootDir, "Machine%i"%MachineNum)
                if not os.path.isdir(MachinePath):
                    os.mkdir(MachinePath)
                    print "creating:", MachinePath
                cfile.CommandFileName = os.path.join(MachinePath, "command.txt")
            RunNum += 1
            #print "Run number is:", RunNum
            start_time = time.rstrip()
            filename = "time" + string.zfill(`i+1`,3)+".nc"
            TrajectoryFile.write(os.path.join(OutputFilePath, filename + '\n'))
            if setup.ReleaseLength == 0:
                end_time = None
            else:
                end_time = batch_gnome.str2DT(start_time) + timedelta(hours=setup.ReleaseLength)
                end_time = batch_gnome.DT2str(end_time)
            Run = batch_gnome.TapRun(start_time,
                                    end_time,
                                    batch_gnome.ConvertToNW(Site),
                                    filename,
                                    OutputFilePath,
                                    Windage=setup.LE_Windage)	
            cfile.AddRun(Run)

TrajectoryFile.close()

if cfile.Runs:
    print "Writing command file:", cfile.CommandFileName
    cfile.write()
        
