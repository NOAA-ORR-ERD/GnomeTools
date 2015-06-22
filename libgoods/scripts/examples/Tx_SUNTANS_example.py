#!/usr/bin/env python
from libgoods import tri_grid, nctools, data_files_dir
import os 
import numpy as np

'''
Example script for Xianlong (TAMU) for working with his SUNTANS model
for Galveston Bay in GNOME

'''

# specify local file or opendap url 
data_file = os.path.join(data_files_dir,'GalvCoarse_2014_0000.nc')

# the utools class requires a mapping of specific model variable names (values)
# to common names (keys) so that the class methods can work with FVCOM, SELFE,
# and ADCIRC which have different variable names
# (This seemed easier than finding them by CF long_names etc)
#!!!!!!!!txsuntans output on server does not include eles_surrounding_ele info
#I have it saved as a netcdf file included in libgoods data_files directory
var_map = { 
            'time':'time',\
            'u_velocity':'uc', \
            'v_velocity':'vc', \
            'nodes_surrounding_ele':'cells',\
            'eles_surrounding_ele':'nbe',\
            'edge_node_connectivity':'edges',\
          }  

# class instantiation creates a netCDF Dataset object as an attribute
txsuntans = tri_grid.ugrid(data_file)

# get longitude, latitude, and time variables
print 'Downloading data dimensions'
txsuntans.get_dimensions(var_map)

# UTM coordinates -- calculate lat/lon
x = txsuntans.Dataset.variables['xp'][:]
y = txsuntans.Dataset.variables['yp'][:]
lon = np.ones_like(x); lat = np.ones_like(x)
for ii in range(len(x)):
    lat[ii], lon[ii] = nctools.utmToLatLng(14,x[ii],y[ii])
txsuntans.data['lon'] = lon
txsuntans.data['lat'] = lat
txsuntans.atts['lon'] = {'long_name': 'longitude'}
txsuntans.atts['lat'] = {'long_name': 'latitude'}

# get grid topo variables (nbe, nv)
print 'Downloading grid topo variables'
txsuntans.get_grid_topo(var_map)

edge_types = txsuntans.Dataset.variables['mark'][:].tolist()
#"0 - computational; 1 - closed; 2 flux BC; 3 - stage BC; 4 - other BC; 5 - interproc; 6 - ghost
#Based on personal communication we use 1,2, and 3
bound_id, bound_type = zip(*[(i,x) for i,x in enumerate(edge_types) if x>0 and x<4])
bound_segs = txsuntans.data['edges'][bound_id,:] + 1 #Node numbering starts at 0, GNOME expects it to start at 1 -- adjust nv and bnd
bound_type_gnome = (np.array(bound_type)<=2).choose(1,0)

txsuntans.order_boundary(bound_segs.tolist(),list(bound_type_gnome))
txsuntans.atts['bnd'] = {'long_name':'Boundary segment information required for GNOME model'} 

#Node numbering starts at 0, GNOME expects it to start at 1 -- adjust nv and bnd
txsuntans.data['nv'] = txsuntans.data['nv'] + 1


## GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
txsuntans.atts['nbe']['order'] = 'ccw'

'''
!!!!!!!!!!!!!All the stuff above here only has to be done once -- if you want to 
process multiple files, I'd put a loop here and just keep overwriting 
txsuntans.data['u'] and ['v'] and incrementing the output file name
Also need to change txsuntans.data['time'] appropriately
'''
# get the data
print 'Loading u/v'
txsuntans.get_data(var_map,zindex=0) 

  
print 'Writing to GNOME file'
txsuntans.write_unstruc_grid(os.path.join(data_files_dir, 'txsuntans_example.nc'))

txsuntans.Dataset.close()