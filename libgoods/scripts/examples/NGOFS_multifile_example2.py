#!/usr/bin/env python
from libgoods import utools, nctools, noaa_coops, data_files_dir
import datetime
import os 
import numpy as np

'''
Sample script to retrieve data from unstructured grid netcdf "file" (can be
OPeNDAP url), generate necessary grid topology (boundary info), and write 
GNOME compatible output.

The boundary file is saved to the data files directory so it only needs 
to be generated once (unless you are subsetting the grid).

To process multiple files (urls) we use
a file list loop after the grid topo vars are loaded (as
this only has to be done once). 
'''

# specify local file or opendap url
sdate = datetime.date(2014,3,21)
edate = datetime.date(2014,3,30)
flist = noaa_coops.make_forecast_list('nwgofs',sdate,edate)

# the utools class requires a mapping of specific model variable names (values)
# to common names (keys) so that the class methods can work with FVCOM, SELFE,
# and ADCIRC which have different variable names
# (This seemed easier than finding them by CF long_names etc)
var_map = { 'longitude':'lon', \
            'latitude':'lat', \
            'time':'time', \
            'u_velocity':'u', \
            'v_velocity':'v', \
            'nodes_surrounding_ele':'nv',\
            'eles_surrounding_ele':'nbe',\
          }  

# class instantiation creates a netCDF Dataset object as an attribute -- 
# use the first file in the list only
nwgofs = utools.ugrid(flist[0])

# get longitude, latitude, and time variables
print 'Downloading data dimensions'
nwgofs.get_dimensions(var_map)

#display available time range for model output
nctools.show_tbounds(nwgofs.Dataset.variables['time'])

# get grid topo variables (nbe, nv)
print 'Downloading grid topo variables'
nwgofs.get_grid_topo(var_map)

# GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
nwgofs.atts['nbe']['order'] = 'cw'


 # GNOME requires boundary info -- this file can be read form data_files directory
# if saved or generated
print 'Loading/generating boundary segments'
ss_bndry_file = os.path.join(data_files_dir, 'nwgofs.bry')
nwgofs.write_bndry_file('nwgofs',ss_bndry_file)
nwgofs.read_bndry_file(ss_bndry_file)
    
# get the data
nwgofs.get_data(var_map)
print 'Downloading data'
for file in flist[1:]:
    print file
    try:
        nwgofs_next = utools.ugrid(file)
        # check time units are consistent among files
        t = nwgofs_next.Dataset.variables[var_map['time']]
        if t.units != nwgofs.atts['time']['units']:
            print 'different time units among files'
            break
        nwgofs_next.data['time'] = t[:]
        nwgofs_next.get_data(var_map)
        # add data to original nwgofs object
        nwgofs.data['time'] = np.concatenate((nwgofs.data['time'],nwgofs_next.data['time']))
        nwgofs.data['u'] = np.concatenate((nwgofs.data['u'],nwgofs_next.data['u']))
        nwgofs.data['v'] = np.concatenate((nwgofs.data['v'],nwgofs_next.data['v']))
    except:
        print 'File not found'
    

    
print 'Writing to GNOME file'
nwgofs.write_unstruc_grid(os.path.join(data_files_dir, 'nwgofs_multifile_example.nc'))