#!/usr/bin/env python
from libgoods import utools, nctools, data_files_dir
import datetime as dt
import os 

# specify local file or opendap url
data_url = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NGOFS/MODELS/201308/nos.ngofs.fields.f000.20130801.t03z.nc'

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
nctools.show_tbounds(ngofs.time)

# get grid topo variables (nbe, nv)
print 'Downloading grid topo variables'
ngofs.get_grid_topo(var_map)

#subsetting
print 'Subsetting'
nl = 31; sl = 26.5
wl = -91; el = -87.5

# Find all nodes and complete elements in subset box, lat/lon subset variables
ngofs.find_nodes_eles_in_ss(nl,sl,wl,el)

# Find subset boundary -- if any of the subset boundary segments correspond to segments
# in the full domain boundary then use boundary type info -- otherwise assume its 
# an open water boundary (then write this new subset boundary to a file)
bndry_file = os.path.join(data_files_dir, 'ngofs.bry') #already exists (no check!)
ngofs.ss_land_bry_segs = ngofs.remap_bry_nodes(bndry_file)
ss_bndry_file = os.path.join(data_files_dir, 'ngofs_ss.bry')
ngofs.write_bndry_file('subset',ss_bndry_file)

# Download u/v -- this is done in multiple OPeNDAP calls of contiguous data blocks
print 'Downloading data'
ngofs.get_data(var_map,nindex=ngofs.nodes_in_ss) #All time steps in file (i.e. tindex=None)

# GNOME requires boundary info -- this file can be read form data_files directory
# if saved or generated
bndry_file = os.path.join(data_files_dir, 'ngofs_ss.bry')
ngofs.read_bndry_file(bndry_file)
 
# GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
ngofs.atts['nbe']['order'] = 'cw'
  
print 'Writing to GNOME file'
ngofs.write_unstruc_grid(os.path.join(data_files_dir, 'ngofs_example_ss.nc'))