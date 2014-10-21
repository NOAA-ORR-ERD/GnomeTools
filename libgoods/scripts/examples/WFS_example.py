#!/usr/bin/env python
from libgoods import utools, nctools, data_files_dir
import datetime as dt
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
data_url = 'http://crow.marine.usf.edu:8080/thredds/dodsC/WFS_FVCOM_NF_model/USF_WFCOM_NF_20140530.nc'

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
wfs = utools.ugrid(data_url)

# get longitude, latitude, and time variables
print 'Downloading data dimensions'
wfs.get_dimensions(var_map)

#display available time range for model output
nctools.show_tbounds(wfs.Dataset.variables['time'])

# get grid topo variables (nbe, nv)
print 'Downloading grid topo variables'
wfs.get_grid_topo(var_map)
wfs.build_face_face_connectivity()
# get the data
print 'Downloading data'
#wfs.get_data(var_map,tindex=[0,1,1]) #First time step only
#wfs.get_data(var_map) #All time steps in file
wfs.get_data(var_map)

# GNOME requires boundary info -- this file can be read form data_files directory
# if saved or generated
bndry_file = os.path.join(data_files_dir, 'WFS.bry')

wfs.write_bndry_file('wfs',bndry_file)
wfs.read_bndry_file(bndry_file)

# GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
wfs.atts['nbe']['order'] = 'ccw'

print 'Writing to GNOME file'
wfs.write_unstruc_grid_only(os.path.join(data_files_dir, 'WFS_grid.nc'))
wfs.write_unstruc_grid(os.path.join(data_files_dir, 'WFS_example.nc'))