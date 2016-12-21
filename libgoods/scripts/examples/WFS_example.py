#!/usr/bin/env python
from __future__ import print_function
from libgoods import tri_grid, nctools, data_files_dir
import os 
import datetime

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
server_stem_url = 'http://crow.marine.usf.edu:8080/thredds/dodsC/WFS_FVCOM_NF_model/'
file_stem = 'USF_WFCOM_NF_' #20140530.nc'
d = datetime.date.today()
data_url = server_stem_url + file_stem + str(d.year) + str(d.month).zfill(2) + str(d.day).zfill(2) + '.nc'

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
          }  

# class instantiation creates a netCDF Dataset object as an attribute
wfs = tri_grid.ugrid(data_url)

# get longitude, latitude, and time variables
print('Downloading data dimensions')
wfs.get_dimensions(var_map)

#display available time range for model output
nctools.show_tbounds(wfs.Dataset.variables['time'])

# get grid topo variables (nbe, nv)
print('Downloading grid topo variables')
wfs.get_grid_topo(var_map)

# find and order the boundary
print('Finding boundary segs')
bnd = wfs.find_bndry_segs()
print('Ordering boundary segs and assigning types')
ow1 = 1; ow2 = 190; #nodes defining start/end of open water boundary
seg_types = []
for b in bnd:
    if max(b) <= ow2 and min(b) >=ow1: #open water
        seg_types.append(1)
    else:
        seg_types.append(0)
wfs.order_boundary(bnd,seg_types)

# get the data
print('Downloading data')
#wfs.get_data(var_map,tindex=[0,1,1]) #First time step only
#wfs.get_data(var_map) #All time steps in file
wfs.get_data(var_map)

# GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
# why ccw?
wfs.atts['nbe']['order'] = 'ccw'

print('Writing to GNOME file')
wfs.write_unstruc_grid_only(os.path.join(data_files_dir, 'WFS_grid.nc'))
wfs.write_unstruc_grid(os.path.join(data_files_dir, 'WFS_example.nc'))
