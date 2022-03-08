#!/usr/bin/env python

"""
script that reads netcdf particle files from GNOME, contours the particles,
makes shape files of the particles, and creates and ERMA data pacakge out
of all of it.

"""
import os
import sys
import shutil
import datetime
from post_gnome import make_layer_file, nc2shape
import yaml

def create_package(params_file):

    params = yaml.safe_load(open(params_file, encoding="utf-8"))

    # make directory structure
    try:
        print("removing:", params['package_dir'])
        shutil.rmtree(params['package_dir'])
    except OSError:
        pass

    os.mkdir(params['package_dir'])
    os.mkdir(os.path.join(params['package_dir'], 'attachments'))
    os.mkdir(os.path.join(params['package_dir'], 'layers'))
    os.mkdir(os.path.join(params['package_dir'], 'source_files'))
    os.mkdir(os.path.join(params['package_dir'], 'support_files'))

    # add attachments
    try:
        if len(params['attachments']) > 0:
            for f in params['attachments']:
                attach_filename = os.path.split(f)[-1]
                print(attach_filename)
                shutil.copy(f, os.path.join(params['package_dir'], 'attachments', attach_filename))
                params['attachment_file'] = attach_filename #need to figure out this for more than one
        else:
            params['attachment_file'] = None
    except:
        params['attachment_file'] = None

    #determine time to convert
    if params['t2convert'] is not None:
        params['t2convert'] = datetime.datetime.strptime(params['t2convert'],'%Y-%m-%d %H:%M')

    # make folder name and append
    if params['folder_name'] is None:
        params['folder_name'] = "Trajectory for " + params['t2convert'].strftime('%Y-%m-%d %H:%M')
    params['folder_path'].append(params['folder_name'])


    if params['styling'] == 'points_simple':
        params['classitem'] = 'status'

        #GNOME best estimate points

        # make shapefiles
        fn = os.path.join(params['gnome_dir'], params['particle_file'])
        traj_zipfname = nc2shape.points(fn, params['package_dir'], params['t2convert'])

        #make layer files
        params['shape_zipfilename'] =  os.path.split(traj_zipfname)[-1]
        print(traj_zipfname)
        try:
            params['title'] = 'Best estimate particles '  + params['t2convert'].strftime('%b %d %Y %H:%M')
        except AttributeError:
            params['title'] = 'Best estimate particles'

        params['color'] = 'black'
        params['color_beached'] = 'black'
        make_layer_file.particles(params['package_dir'],'traj',params)

        # Red uncertainty points if params['uncertain']=True
        if params['uncertain']:

            #make shapefiles
            ufn = os.path.join(params['gnome_dir'],params['particle_file'].split('.')[0] + '_uncertain.nc')
            uncert_zipfname = nc2shape.points(ufn,params['package_dir'],params['t2convert'])

            #make layer file
            params['shape_zipfilename'] =  os.path.split(uncert_zipfname)[-1]
            params['color'] = 'red'
            params['color_beached'] = 'red'
            try:
                params['title'] = 'Uncertainty particles '  + params['t2convert'].strftime('%b %d %Y %H:%M')
            except AttributeError:
                params['title'] = 'Uncertainty particles'
            params['attachment_file'] = None
            make_layer_file.particles(params['package_dir'],'uncert',params)

    elif params['styling'] == 'points_forecast':
        params['classitem'] = 'surf_conc'

        # ***************Best estimate particles

         # make shapefiles
        fn = os.path.join(params['gnome_dir'], params['particle_file'])
        traj_zipfname = nc2shape.points(fn,params['package_dir'],params['t2convert'],status_code=2)

        #make layer files
        params['shape_zipfilename'] =  os.path.split(traj_zipfname)[-1]
        try:
            params['title'] = 'Floating Oil at '  + params['t2convert'].strftime('%b %d %Y %H:%M')
        except AttributeError:
            params['title'] = 'Floating Oil '
        make_layer_file.particles(params['package_dir'],'traj',params)

       # ***************Beached particles

        # make shapefiles
        fn = os.path.join(params['gnome_dir'], params['particle_file'])
        traj_zipfname = nc2shape.points(fn,params['package_dir'],params['t2convert'],status_code=3, shapefile_name=params['particle_file'].split('.nc')[0]+'_beached')

        #make layer files
        params['shape_zipfilename'] =  os.path.split(traj_zipfname)[-1]
        try:
            params['title'] = 'Beached particles at '  + params['t2convert'].strftime('%b %d %Y %H:%M')
        except AttributeError:
            params['title'] = 'Beached particles '
        params['color'] = 'red'
        params['color_beached'] = 'red'
        params['classitem'] = 'status'
        make_layer_file.particles(params['package_dir'],'beached',params)
       # **************Uncertainty contour if params['uncertain']=True

        if params['uncertain']:
            # make shapefile
            ufn = os.path.join(params['gnome_dir'], params['particle_file'].split('.')[0] + '_uncertain.nc')
            uncert_zipfname = nc2shape.contours(ufn,
                                                params['package_dir'],
                                                params['t2convert'],
                                                levels=[0.01],
                                                names=['Uncertainty'],
                                                include_beached=True
                                                )
            # make layer file
            params['shape_zipfilename'] =  os.path.split(uncert_zipfname)[-1]
            try:
                params['title'] = 'Uncertainty contour ' + params['t2convert'].strftime('%b %d %Y %H:%M')
            except AttributeError:
                params['title'] = 'Uncertainty contour'
            params['attachment_file'] = None
            make_layer_file.contours(params['package_dir'],'uncert',params)

    elif params['styling'] == 'contours_forecast':

        # Best estimate contours

        # make shapefiles
        fn = os.path.join(params['gnome_dir'], params['particle_file'])
        traj_zipfname = nc2shape.contours(fn, params['package_dir'], params['t2convert'])
        # make layer files
        params['shape_zipfilename'] =  os.path.split(traj_zipfname)[-1]
        try:
            params['title'] = 'Best estimate contours ' + params['t2convert'].strftime('%b %d %Y %H:%M')
        except AttributeError:
            params['title'] = 'Best estimate contours'
        params['color'] = 'black'
        #params['SinglePoly'] = True
        make_layer_file.contours(params['package_dir'],'traj',params)

        if params['uncertain']:
            # make shapefile
            ufn = os.path.join(params['gnome_dir'], params['particle_file'].split('.')[0] + '_uncertain.nc')
            uncert_zipfname = nc2shape.contours(ufn,
                                                params['package_dir'],
                                                params['t2convert'],
                                                levels=[0.1],
                                                names=['Uncertainty']
                                                )
            # make layer file
            params['shape_zipfilename'] =  os.path.split(uncert_zipfname)[-1]
            try:
                params['title'] = 'Uncertainty contour ' + params['t2convert'].strftime('%b %d %Y %H:%M')
            except AttributeError:
                params['title'] = 'Uncertainty contour'
            params['attachment_file'] = None
            make_layer_file.contours(params['package_dir'],'uncert',params)

        # Add beached particles
        params['classitem'] = 'status'
        # make shapefiles
        traj_zipfname = nc2shape.points(fn, params['package_dir'], params['t2convert'], beached_only=True, shapefile_name = 'Beached')
        #make layer files
        params['shape_zipfilename'] =  os.path.split(traj_zipfname)[-1]
        try:
            params['title'] = 'Beached particles '  + params['t2convert'].strftime('%b %d %Y %H:%M')
        except AttributeError:
            params['title'] = 'Beached particles'
        params['color'] = 'black'
        params['color_beached'] = 'black'
        make_layer_file.particles(params['package_dir'],'beached',params)

    else:
        print('Must specify either points or contours')

    shutil.make_archive(params['package_dir'],'zip',root_dir=params['package_dir'])

USAGE = """
make_erma_data_package.py params.yml

params.yml holds the configuration for how you want the data packet configured

if you don't pass in a config file -- this script will dump a sample one to update
"""

EXAMPLE_PARAMS_FILE = """
gnome_dir : .
particle_file : Guam_particles2.nc
uncertain : True # include uncertainty contour (boolean)

# If not null -- the one timestep to output
t2convert : 2019-03-12 03:00 # null for animation time series or %Y-%m-%d %H:%M
metadata : [Hypothetical trajectory for Guam Science of Oil Spills Training Class, March 11-15, 2019]

# you can attach arbitrary other files -- PDFs, etc.
attachments : ['file1','file2'] # relative path 2 single file OR list of paths or []

# This will be a name of the resulting zipfile
package_dir : 24hr_fc # name of directory to create files in (will be created and destroyed if need be

# Styling options:
#   points_simple (floating/beached)
#   points_forecast (surface_concentation)
#   contours_forecast (only the contours, not points)
styling : points_forecast  # points_simple (floating/beached); points_forecast (surface_concentation); or contours_forecast

# ERMA Site name
site_name : pacific
event : Guam SOS Training Scenario (DRILL)

# ERMA's folder structure
folder_path : [Incidents & Drills, Guam SOS/SCAT Training Class Scenario, Trajectories]
folder_name : 12 hour trajectory # if null folder will be named "Trajectory for %Y-%m-%d %H:%M"
"""


if __name__ == "__main__":
    try:
        params_file = sys.argv[1]
        create_package(params_file)
    except IndexError:
        print(USAGE)
        with open("example_params.yml", 'w', encoding='utf-8') as outfile:
            outfile.write(EXAMPLE_PARAMS_FILE)



