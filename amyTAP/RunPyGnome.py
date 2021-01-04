#!/usr/bin/env python

import os
from datetime import datetime, timedelta
import gc

from gnome.spill import point_line_release_spill
from gnome.outputters import Renderer, NetCDFOutput
from gnome.model import Model
from gnome.map import MapFromBNA
from gnome.movers import PyCurrentMover, GridWindMover, RandomMover, GridCurrentMover, WindMover, CatsMover, constant_wind_mover
from gnome.environment import Wind, Tide

from TAP_Setup import setup

# create command file dir if it doesn't exist
print "RootDir is", setup.RootDir

def make_dir(path):
    #path = os.path.join(setup.RootDir, name)
    if not os.path.exists(path):
        os.makedirs(path)

def setup_model():
    print 'initializing the model'
    # start with default time,duration...this will be changed when model is run
    model = Model() #change to use all defaults and set time_step also in Setup_TAP!!
    mapfile = os.path.join(setup.MapFileDir, setup.MapFileName)
    print 'adding the map: ', mapfile
    model.map = MapFromBNA(mapfile, refloat_halflife=0.0)  # seconds
    
    print 'adding a GridCurrentMover:'
    if setup.curr_fn.endswith('.nc'):
        c_mover = GridCurrentMover(filename=setup.curr_fn,extrapolate=True)
    elif setup.curr_fn.endswith('.cur'):
        tide = Tide(setup.tide_fn)
        c_mover = CatsMover(filename=setup.curr_fn,tide=tide)
    model.movers += c_mover

    print 'adding a WindMover:'
    if setup.wind_fn is not None:
        w = Wind(filename=setup.wind_fn)
        w_mover = WindMover(w)
    elif setup.wind_data is not None:
        w_mover = constant_wind_mover(setup.wind_data[0],setup.wind_data[1],units='knots')    
    # w_mover = GridWindMover(wind_file=setup.w_filelist)
    model.movers += w_mover

    if setup.diff_coef is not None:
        print 'adding a RandomMover:'
        random_mover = RandomMover(diffusion_coef=setup.diff_coef) #in cm/s 
        model.movers += random_mover

    return model


# load up the start positions
start_positions = open(os.path.join(setup.RootDir,
                                    setup.CubeStartSitesFilename)).readlines()
start_positions = [pos.split(',') for pos in start_positions]
start_positions = [( float(pos[0]), float(pos[1]) ) for pos in start_positions]

# model timing
release_duration = timedelta(hours=setup.ReleaseLength)
run_time = timedelta(hours=setup.TrajectoryRunLength)

model = setup_model()
model.duration = run_time
    
# # loop through seasons
for Season in setup.StartTimeFiles:
    SeasonName = Season[1]
    start_times = open(Season[0],'r').readlines()[:setup.NumStarts]

    make_dir(os.path.join(setup.RootDir,setup.TrajectoriesPath,SeasonName))
    print Season

    # get and parse start times in this season
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

    # loop through start times
    for time_idx in range(0,setup.NumStarts): 
        gc.collect()    # garbage collection

        start_time = start_dt[time_idx]
        end_time = start_time + run_time
        print start_time, end_time  

        for pos_idx, start_position in enumerate(start_positions):
        # for pos_idx in setup.RunSites:
            start_position = start_positions[pos_idx]

            OutDir = os.path.join(setup.RootDir,setup.TrajectoriesPath,SeasonName,'pos_%03i'%(pos_idx+1))
            make_dir(OutDir)

            print pos_idx, time_idx
            print "Running: start time:", start_time,
            print "At start location:",   start_position

            ## set the spill to the location
            spill = point_line_release_spill(num_elements=setup.NumLEs,
                                             start_position=( start_position[0], start_position[1], 0.0 ),
                                             release_time=start_time,
                                             end_release_time=start_time+release_duration
                                             )

            # set up the renderer
            image_dir = os.path.join(setup.ImagesPath, SeasonName, 'images_pos_%03i-time_%03i'%(pos_idx+1, time_idx))
            renderer = Renderer(os.path.join(setup.MapFileDir, setup.MapFileName),
                                image_dir,
                                image_size=(800, 600),
                                output_timestep=timedelta(hours=6))
            make_dir(image_dir)

            # setup netcdf
            netcdf_output_file = os.path.join(OutDir,
                                              'pos_%03i-t%03i_%08i.nc'%(pos_idx+1, time_idx,
                                                int(start_time.strftime('%y%m%d%H'))),
                                              )


            model.start_time = start_time

            ## clear the old outputters
            model.outputters.clear()
            model.outputters += renderer
            model.outputters += NetCDFOutput(netcdf_output_file,output_timestep=timedelta(hours=setup.GnomeOutputTimestepHours))

            # clear out the old spills:
            model.spills.clear()
            model.spills+=spill

            model.full_run(rewind=True)
 

