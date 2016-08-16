#!/usr/bin/env python

import os
import shutil
import datetime
import nc2shape
import make_layer_file
reload(make_layer_file)

# GNOME file info***************************************************************
gnome_dir = ''
metadata = """Data sources:

Overflights from 5/16/2016
No recoverable oil observed.

Currents:
US Navy American Seas Model

Winds:
NWS GFS
"""

uncertain = True
outfile = 'Forecast-5-16.nc'
t2convert = datetime.datetime(2016, 5, 19, 8, 0)
# attachments = os.path.join(gnome_dir,'images','ChesapeakeBay_anim.gif')
# Name for package directory****************************************************
package_dir = 'GreenCanyon_ERMA_forecast_5-16_48hr'
# ERMA layer info
plot_type = 'contours'  # points or contours
params = {}
params['site_name'] = "gulfofmexico"
params['event'] = "Green Canyon 248 12-May-2015"
params['metadata'] = metadata
params['folder_path'] = ['Green Canyon 248 Incident',
                         'Trajectory Forecast',
                         'Trajectory for ' + t2convert.strftime('%b %d %Y %H:%M')]
# ******************************************************************************

# make directory structure
try:
    print "removing:", package_dir
    shutil.rmtree(package_dir)
except OSError:
    pass

os.mkdir(package_dir)
os.mkdir(os.path.join(package_dir, 'attachments'))
os.mkdir(os.path.join(package_dir, 'layers'))
os.mkdir(os.path.join(package_dir, 'source_files'))
os.mkdir(os.path.join(package_dir, 'support_files'))

# add attachments
try:
    attach_file = os.path.split(attachments)[-1]
    shutil.copy(attachments, os.path.join(package_dir, 'attachments', attach_file))
    params['attachment_file'] = attach_file
except:
    print 'No attachments'
    params['attachment_file'] = None

if plot_type == 'points':
    # make shapefiles
    fn = os.path.join(gnome_dir, outfile)
    traj_zipfname = nc2shape.points(fn, package_dir,t2convert)
    if uncertain:
        ufn = os.path.join(gnome_dir,outfile.split('.')[0] + '_uncertain.nc')
        uncert_zipfname = nc2shape.points(ufn,package_dir,t2convert)

    #make layer files
    params['shape_zipfilename'] = traj_zipfname
    params['title'] = 'Best estimate particles: ' + t2convert.strftime('%b %d %Y %H:%M')
    params['color'] = 'black'
    make_layer_file.particles(package_dir,'traj',params)

    if uncertain:
        params['shape_zipfilename'] = uncert_zipfname
        params['color'] = 'red'
        params['title'] = 'Uncertainty particles: ' + t2convert.strftime('%b %d %Y %H:%M')
        params['attachment_file'] = None
        make_layer_file.points(package_dir,'uncert',params)

elif plot_type == 'contours':
    #make shapefiles
    fn = os.path.join(gnome_dir,outfile)
    traj_zipfname = nc2shape.contours(fn,package_dir,t2convert)
    if uncertain:
        ufn = os.path.join(gnome_dir,outfile.split('.')[0] + '_uncert.nc')
        print ufn
        uncert_zipfname = nc2shape.contours(ufn,package_dir,t2convert,levels=[0.1,],names=['Uncertainty',])
        
    #make layer files
    params['shape_zipfilename'] = traj_zipfname
    params['title'] = 'Best estimate contours: ' + t2convert.strftime('%b %d %Y %H:%M')
    params['color'] = 'black'
    params['SinglePoly'] = True
    make_layer_file.contours(package_dir,'traj',params)
    
    if uncertain:
        params['shape_zipfilename'] = uncert_zipfname
        params['title'] = 'Uncertainty contour: ' + t2convert.strftime('%b %d %Y %H:%M')
        params['attachment_file'] = None
        params['Uncertain'] = True
        make_layer_file.contours(package_dir,'uncert',params)

else:
    print 'Must specify either points or contours'
    
shutil.make_archive(package_dir,'zip',root_dir=package_dir)
