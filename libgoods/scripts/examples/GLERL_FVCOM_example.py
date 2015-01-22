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
data_url = 'http://tds.glos.us/thredds/dodsC/FVCOM/SLRFVM-Latest-Forecast.nc'
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
slrfvm = tri_grid.ugrid(data_url)

# get longitude, latitude, and time variables
print 'Downloading data dimensions'
slrfvm.get_dimensions(var_map)

#display available time range for model output
nctools.show_tbounds(slrfvm.Dataset.variables['time'])

# get grid topo variables (nbe, nv)
print 'Downloading grid topo variables'
slrfvm.get_grid_topo(var_map)
# GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
slrfvm.atts['nbe']['order'] = 'cw'

# find and order the boundary
print 'Finding boundary'
bnd = slrfvm.find_bndry_segs()
print 'Ordering boundary'
seg_types = [0] * len(bnd)
slrfvm.order_boundary(bnd,seg_types)

# get the data
print 'Downloading data'
#slrfvm.get_data(var_map,tindex=[0,1,1]) #First time step only
slrfvm.get_data(var_map) #All time steps in file
 
print 'Writing to GNOME file'
#slrfvm.write_unstruc_grid_only(os.path.join(data_files_dir, 'slrfvm_grid.nc'))
slrfvm.write_unstruc_grid(os.path.join(data_files_dir, 'slrfvm_example.nc'))