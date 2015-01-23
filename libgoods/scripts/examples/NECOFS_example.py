#!/usr/bin/env python
from libgoods import tri_grid, nctools, data_files_dir
reload(tri_grid)
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
#data_url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/archives/necofs_mb'
data_url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc'
grid = 'massb' #gom2, gom3, or massb
#data_url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc' 
# the tri_grid class requires a mapping of specific model variable names (values)
# to common names (keys) so that the class methods can work with FVCOM, SELFE,
# and ADCIRC which have different variable names
# (This seemed easier than finding them by CF long_names etc)
var_map = { 'longitude':'lon', \
            'latitude':'lat', \
            'time':'time', \
            'u_velocity':'ua', \
            'v_velocity':'va', \
            'nodes_surrounding_ele':'nv',\
            'eles_surrounding_ele':'nbe',\
          }  

# class instantiation creates a netCDF Dataset object as an attribute
necofs = tri_grid.ugrid(data_url)

# get longitude, latitude, and time variables
print 'Downloading data dimensions'
necofs.get_dimensions(var_map)

#display available time range for model output
nctools.show_tbounds(necofs.Dataset.variables['time'])

# determine subset indices (temporal)
print 'Determining subset indices'
# start = dt.datetime(2011,4,13,6,0,0)
# stop = dt.datetime(2011,4,13,6,0,0)
# tindex = nctools.get_tindex(necofs.time,start,stop)

# get grid topo variables (nbe, nv)
print 'Downloading grid topo variables'
necofs.get_grid_topo(var_map)

# find and order the boundary
print 'Finding boundary segs'
bnd = necofs.find_bndry_segs()
print 'Ordering boundary segs and assigning types'
if grid.lower() == 'gom2':
    ow1 = 1; ow2 = 60;
elif grid.lower() == 'gom3':
    ow1 = 1; ow2 = 120;
elif grid.lower() == 'massb':
    ow1 = 1; ow2 = 124;
seg_types = []
for b in bnd:
    if max(b) <= ow2 and min(b) >=ow1: #open water
        seg_types.append(1)
    else:
        seg_types.append(0)
necofs.order_boundary(bnd,seg_types)

## get the data
#print 'Downloading data'
#necofs.get_data(var_map,tindex=[0,1,1]) #First time step only
necofs.get_data(var_map) #All time steps in file
#necofs.get_data(var_map,tindex=tindex)
#

## GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
necofs.atts['nbe']['order'] = 'cw'
#
#print 'Writing to GNOME file'
necofs.write_unstruc_grid(os.path.join(data_files_dir, 'NECOFS_massb_example.nc'))