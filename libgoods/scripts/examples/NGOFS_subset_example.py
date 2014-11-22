#!/usr/bin/env python
from libgoods import utools, nctools, data_files_dir
import os 

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
#data_url = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NGOFS/MODELS/201403/nos.ngofs.fields.f000.20140324.t09z.nc'
data_url = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NGOFS/MODELS/201405/nos.ngofs.fields.n000.20140501.t09z.nc'
#data_url = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NGOFS/MODELS/201403/nos.ngofs.fields.nowcast.20130301.t03z.nc'

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

# class instantiation creates a netCDF Dataset object as an attribute
ngofs = utools.ugrid(data_url)

# get longitude, latitude, and time variables
print 'Downloading data dimensions'
ngofs.get_dimensions(var_map)

#display available time range for model output
nctools.show_tbounds(ngofs.Dataset.variables['time'])

# get grid topo variables (nbe, nv)
print 'Downloading grid topo variables'
ngofs.get_grid_topo(var_map)
# GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
ngofs.atts['nbe']['order'] = 'cw'

#subsetting
print 'Subsetting'
nl = 29.7; sl = 28.1
wl = -97; el = -94


# Find all nodes and complete elements in subset box, lat/lon subset variables
ngofs.find_nodes_eles_in_ss(nl,sl,wl,el)
print ngofs.nodes_in_ss
print ngofs.eles_in_ss

# GNOME requires boundary info -- this file can be read form data_files directory
# if already generated (!!!for this particular subset!!!)
# Find subset boundary -- if any of the subset boundary segments correspond to segments
# in the full domain boundary then use boundary type info -- otherwise assume its 
# an open water boundary (then write this new subset boundary to a file)

bndry_file = os.path.join(data_files_dir, 'ngofs.bry') #already exists (no check!)
ngofs.ss_land_bry_segs = ngofs.remap_bry_nodes(bndry_file)
ngofs.write_bndry_file('subset')


# Download u/v -- this is done in multiple OPeNDAP calls of contiguous data blocks
print 'Downloading data'
ngofs.get_data(var_map,nindex=ngofs.nodes_in_ss) #All time steps in file (i.e. tindex=None)
  
print 'Writing to GNOME file'
ngofs.write_unstruc_grid(os.path.join(data_files_dir, 'ngofs_ss_may.nc'))