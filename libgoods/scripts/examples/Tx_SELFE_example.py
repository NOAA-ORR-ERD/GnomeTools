#!/usr/bin/env python
from libgoods import tri_grid, nctools, data_files_dir
import datetime as dt
import os 
import numpy as np
from netCDF4 import date2num

'''
Sample script to retrieve data from unstructured grid netcdf "file" (can be
OPeNDAP url), generate necessary grid topology (boundary info), and write 
GNOME compatible output.

The boundary file is saved to the data files directory so it only needs 
to be generated once (unless you are subsetting the grid).

To process multiple files (urls) either
a) pass the filenames/urls in as a list -- this creates a netcdf4 MFDataset and is
a good option for not too many files (all output is written to one nc file for GNOME 
in this case)
b) add a file list loop -- in this case put it after the grid topo vars are loaded (as
this only has to be done once). See NGOFS_multifile_example.py

'''

# specify local file or opendap url 
data_file = os.path.join(data_files_dir,'tx_selfe_native_output.nc')

# the utools class requires a mapping of specific model variable names (values)
# to common names (keys) so that the class methods can work with FVCOM, SELFE,
# and ADCIRC which have different variable names
# (This seemed easier than finding them by CF long_names etc)
#!!!!!!!!txselfe output on server does not include eles_surrounding_ele info
#I have it saved as a netcdf file included in libgoods data_files directory
var_map = { 'longitude':'lon', \
            'latitude':'lat', \
            'time':'', \
            'u_velocity':'u', \
            'v_velocity':'v', \
            'nodes_surrounding_ele':'ele',\
            'eles_surrounding_ele':'',\
          }  

# class instantiation creates a netCDF Dataset object as an attribute
txselfe = tri_grid.ugrid(data_file)

# get longitude, latitude, and time variables
print 'Downloading data dimensions'
# Normally this would be a call to txself.get_dimensions(var_map) but
# this file has no time info, and only UTM coordinates -- so we do it
# manually

#Enter actual model time here for particular file
model_time = dt.datetime(2013,8,21,0,0,0) #I made this time up
t_units = 'hours since 2012-01-01 00:00:00'
txselfe.data['time'] = [date2num(model_time,t_units),]
txselfe.atts['time'] = {'units':t_units}

x = txselfe.Dataset.variables['x'][:]
y = txselfe.Dataset.variables['y'][:]
lon = np.ones_like(x); lat = np.ones_like(x)
for ii in range(len(x)):
    lat[ii], lon[ii] = nctools.utmToLatLng(14,x[ii],y[ii])
txselfe.data['lon'] = lon
txselfe.data['lat'] = lat
txselfe.atts['lon'] = {'long_name': 'longitude'}
txselfe.atts['lat'] = {'long_name': 'latitude'}

# get grid topo variables (nbe, nv)
print 'Downloading grid topo variables'
txselfe.get_grid_topo(var_map)


# find and order the boundary
print 'Finding boundary segs'
bnd = txselfe.find_bndry_segs()
print 'Ordering boundary segs and assigning types'
ow1 = 1; ow2 = max(max(bnd)); #nodes defining start/end of open water boundary
seg_types = []
for b in bnd:
    if max(b) <= ow2 and min(b) >=ow1: #open water
        seg_types.append(0)
    else:
        seg_types.append(1)
txselfe.order_boundary(bnd,seg_types)

# GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
txselfe.atts['nbe']['order'] = 'ccw'

'''
!!!!!!!!!!!!!All the stuff above here only has to be done once -- if you want to 
process multiple files, I'd put a loop here and just keep overwriting 
txselfe.data['u'] and ['v'] and incrementing the output file name
Also need to change txselfe.data['time'] appropriately
'''
# get the data
print 'Loading u/v'
# txselfe.get_data(var_map,zindex=-1) #First time step only
# make up some data since there is none in the nc file
num_nodes = len(txselfe.data['lon'])
txselfe.data['u'] = np.ones([1,num_nodes],)
txselfe.data['v'] = np.ones([1,num_nodes],)
txselfe.atts['u'] = {'long_name':'eastward_velocity','units':'m/s'}
txselfe.atts['v'] = {'long_name':'northward_velocity','units':'m/s'}
  
print 'Writing to GNOME file'
txselfe.write_unstruc_grid(os.path.join(data_files_dir, 'txselfe_example.nc'))