#!/usr/bin/env python
from libgoods import tri_grid, noaa_coops, nctools, data_files_dir
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
this only has to be done once). See COOPS_FVCOM_multifile_example.py

'''


# specify local file or opendap url -- in this case files are one time step, not aggregated
data_url = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/CREOFS/MODELS/201412/nos.creofs.fields.n000.20141225.t03z.nc'

# the utools class requires a mapping of specific model variable names (values)
# to common names (keys) so that the class methods can work with FVCOM, SELFE,
# and ADCIRC which have different variable names
# (This seemed easier than finding them by CF long_names etc)
#!!!!!!!!CREOFS output on server does not include eles_surrounding_ele info
#I have it saved as a netcdf file included in libgoods data_files directory
var_map = { 'longitude':'lon', \
            'latitude':'lat', \
            'time':'time', \
            'u_velocity':'u', \
            'v_velocity':'v', \
            'nodes_surrounding_ele':'ele',\
            'eles_surrounding_ele':'',\
          }  

# class instantiation creates a netCDF Dataset object as an attribute
# creofs = utools.ugrid(flist) #multiple files
creofs = tri_grid.ugrid(data_url) #single file output

# get longitude, latitude, and time variables
print 'Downloading data dimensions'
creofs.get_dimensions(var_map)

#display available time range for model output
nctools.show_tbounds(creofs.Dataset.variables['time'])

# get grid topo variables (nbe, nv)
print 'Downloading grid topo variables'
creofs.get_grid_topo(var_map)
creofs.atts['nbe']['order'] = 'ccw'

# find and order the boundary
print 'Finding boundary'
bnd = creofs.find_bndry_segs()
print 'Ordering boundary'
seg_types = noaa_coops.specify_bnd_types('creofs',bnd)
creofs.order_boundary(bnd,seg_types)

# get the data
print 'Downloading data'
#creofs.get_data(var_map,tindex=[0,1,1]) #First time step only
creofs.get_data(var_map,zindex=-1) #All time steps in file
 
print 'Writing to GNOME file'
creofs.write_unstruc_grid(os.path.join(data_files_dir, 'creofs_example.nc'))