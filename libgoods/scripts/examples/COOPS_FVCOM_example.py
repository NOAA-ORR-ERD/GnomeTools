#!/usr/bin/env python
from libgoods import tri_grid, nctools, noaa_coops, data_files_dir
import os 

'''
Sample script to retrieve data from NOAA CO-OPS FVCOM netcdf "file" (can be
OPeNDAP url), generate necessary grid topology (boundary info), and write 
GNOME compatible output.

'''
# specify local file or opendap url
data_url = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NWGOFS/MODELS/201501/nos.nwgofs.fields.n000.20150115.t09z.nc'
#data_url = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NWGOFS/MODELS/201501/nos.nwgofs.fields.f001.20150102.t15z.nc'
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
nwgofs = tri_grid.ugrid(data_url)

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

# find and order the boundary
print 'Finding boundary'
bnd = nwgofs.find_bndry_segs()
print 'Ordering boundary'
seg_types = noaa_coops.specify_bnd_types('nwgofs',bnd)
nwgofs.order_boundary(bnd,seg_types)

# get the data
print 'Downloading data'
#nwgofs.get_data(var_map,tindex=[0,1,1]) #First time step only
nwgofs.get_data(var_map) #All time steps in file
 
print 'Writing to GNOME file'
nwgofs.write_unstruc_grid(os.path.join(data_files_dir, 'nwgofs_example.nc'))