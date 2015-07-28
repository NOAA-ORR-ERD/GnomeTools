#!/usr/bin/env python

import os, sys, string, math
from datetime import datetime, timedelta
import importlib

import gnome
from gnome.spill import point_line_release_spill
from gnome.outputters import Renderer, NetCDFOutput

from TAP_Setup import setup

#from batch_gnome import batch_gnome

# create command file dir if it doesn't exist
print "RootDir is", setup.RootDir


# print "Setting up Master Command File"
# cfile = batch_gnome.CommandFile()
# cfile.SaveFilePath   = os.curdir
# cfile.SaveFileName   = setup.SaveFileName
# cfile.RunLength      = setup.TrajectoryRunLength
# cfile.NumLEs         = setup.NumLEs
# cfile.ModelTimeStep  = setup.ModelTimeStep

# NumCubes = len(setup.CubeStartSites)
NumStarts = setup.NumStarts
NumSeasons = len(setup.StartTimeFiles)

## load up our run script
sys.path.append(setup.RootDir)
script = importlib.import_module(setup.PyGnome_script)

print script
print dir(script)

model = script.make_model()

start_positions = open(os.path.join(setup.RootDir,
                                    setup.CubeStartSitesFilename)).readlines()
start_positions = [pos.split(',') for pos in start_positions]
start_positions = [( float(pos[0]), float(pos[1]) ) for pos in start_positions]

start_times = open(os.path.join(setup.RootDir, "All_yearStarts.txt")).readlines()

start_dt = []
for start_time in start_times:
    start_time = [int(i) for i in start_time.split(',')]
    start_time = datetime(start_time[2],
                          start_time[1],
                          start_time[0],
                          start_time[3],
                          start_time[4],
                          )
    start_dt.append(start_time)

run_time = timedelta(hours=setup.TrajectoryRunLength)
model.duration = run_time

# spill = point_line_release_spill(num_elements=setup.NumLEs,
#                                  start_position=(0.0, 0.0, 0.0),
#                                  release_time=start_dt[0],
#                                  end_release_time=start_time+run_time
#                                  )

# model.spills += spill
# model.outputters

def make_dir(name):
    path = os.path.join(setup.RootDir, name)
    if not os.path.exists(path):
        os.makedirs(path)

make_dir(setup.TrajectoriesPath)

for pos_idx, start_position in enumerate(start_positions):
    for time_idx, start_time in enumerate(start_dt):
        print "Running: start time:", start_time,
        print "At start location:",   start_position

        ## set the spill to the location
        spill = point_line_release_spill(num_elements=setup.NumLEs,
                                         start_position=( start_position[0], start_position[1], 0.0 ),
                                         release_time=start_time,
                                         end_release_time=start_time+run_time
                                         )

        # set up the renderer location
        image_dir = os.path.join(setup.RootDir, 'images_pos_%03i-time_%03i.nc'%(pos_idx, time_idx))
        renderer = Renderer(os.path.join(setup.RootDir, setup.MapFileName),
                            image_dir,
                            image_size=(800, 600),
                            output_timestep=timedelta(hours=4))
        #renderer.viewport = ((-120.6666, 33.75),(-119.25, 34.5)) 
        make_dir(image_dir)

        # print "adding netcdf output"
        netcdf_output_file = os.path.join(setup.RootDir,
                                          setup.TrajectoriesPath,
                                          'pos_%03i-time_%03i.nc'%(pos_idx, time_idx),
                                          )


        model.start_time = start_time

        ## clear the old outputters
        model.outputters.clear()
        model.outputters += renderer
        model.outputters += NetCDFOutput(netcdf_output_file)

        # clear out the old spills:
        model.spills.clear()
        model.spills+=spill

        model.full_run(rewind=True)
        # for i, step in enumerate(model):
        #     print i, step
        #     print
        #     for sc in model.spills.items():
        #         print "status_codes:", sc['status_codes']
        #         print "positions:", sc['positions']
        #         print "lw positions:", sc['last_water_positions']




## Old code for command file processing:


# RunsPerMachine = int( math.ceil(float((NumCubes * NumSeasons * NumStarts)) / setup.NumTrajMachines ) )
# print "RunsPerMachine", RunsPerMachine
    


#CubesPerMachine = int( math.ceil(float(len(setup.CubeStartSites))  * NumSeasons / setup.NumTrajMachines) )



# TrajectoryFile = file(os.path.join(setup.RootDir, 'ExpectedTrajectories.txt'),'w')
# MachineNum = 0
# RunNum = 0





# for Season in setup.StartTimeFiles:
#     SeasonName = Season[1]
#     StartTimes = open(Season[0],'r').readlines()[:setup.NumStarts]
#     CubeNum = 0
#     for Site in setup.CubeStartSites:
#         CubeNum += 1
#         OutputFilePath  = os.path.abspath(os.path.join(setup.RootDir, setup.TrajectoriesPath, SeasonName, str(CubeNum).zfill(3) ))
#         #print "OutputFilePath", OutputFilePath
#         #print os.path.abspath(OutputFilePath)
#         #raise Exception("stopping here")
#         for i,time in zip(range(len(StartTimes)),StartTimes):
#             if (RunNum % RunsPerMachine) == 0: # make a new command file
#                 MachineNum += 1
#                 if cfile.Runs:
#                     cfile.write()
#                     cfile.ClearRuns()
#                 MachinePath = os.path.join(setup.RootDir, "Machine%i"%MachineNum)
#                 if not os.path.isdir(MachinePath):
#                     os.mkdir(MachinePath)
#                     print "creating:", MachinePath
#                 cfile.CommandFileName = os.path.join(MachinePath, "command.txt")
#             RunNum += 1
#             #print "Run number is:", RunNum
#             start_time = time.rstrip()
#             filename = "time" + string.zfill(`i+1`,3)+".nc"
#             TrajectoryFile.write(os.path.join(OutputFilePath, filename + '\n'))
#             if setup.ReleaseLength == 0:
#                 end_time = None
#             else:
#                 end_time = batch_gnome.str2DT(start_time) + timedelta(hours=setup.ReleaseLength)
#                 end_time = batch_gnome.DT2str(end_time)
#             Run = batch_gnome.TapRun(start_time,
#                                     end_time,
#                                     batch_gnome.ConvertToNW(Site),
#                                     filename,
#                                     OutputFilePath,
#                                     Windage=setup.LE_Windage)	
#             cfile.AddRun(Run)

# TrajectoryFile.close()

# if cfile.Runs:
#     print "Writing command file:", cfile.CommandFileName
#     cfile.write()
        
