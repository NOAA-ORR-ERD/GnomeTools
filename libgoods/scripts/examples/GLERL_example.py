#!/usr/bin/env python
from libgoods import utools, nctools, data_files_dir
import os 

'''
Sample script to retrieve data from GLERL FVCOM

'''
# specify local file or opendap url
data_file = os.path.join(data_files_dir,'201421600.nc')

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
erie = utools.ugrid(data_file)

# get longitude, latitude, and time variables
print 'Downloading data dimensions'
erie.get_dimensions(var_map)

#display available time range for model output
nctools.show_tbounds(erie.Dataset.variables['time'])

# get grid topo variables (nbe, nv)
print 'Downloading grid topo variables'
erie.get_grid_topo(var_map)

print 'Downloading grid topo variables'
erie.get_grid_topo(var_map)
# GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
erie.atts['nbe']['order'] = 'ccw'

# GNOME requires boundary info -- this file can be read form data_files directory
# if saved or generated
print 'Loading/generating boundary segments'
bndry_file = os.path.join(data_files_dir, 'erie.bry')
try:
    erie.read_bndry_file(bndry_file)
except IOError:
    erie.write_bndry_file('erie',bndry_file)
    erie.read_bndry_file(bndry_file)

# get the data
print 'Downloading data'
#erie.get_data(var_map,tindex=[0,1,1]) #First time step only
erie.get_data(var_map) #All time steps in file
 
print 'Writing to GNOME file'
erie.write_unstruc_grid(os.path.join(data_files_dir, 'erie_example.nc'))