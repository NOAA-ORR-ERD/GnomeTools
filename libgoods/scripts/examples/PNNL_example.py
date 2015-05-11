#!/usr/bin/env python
from libgoods import tri_grid, nctools
reload(tri_grid)
import os 
import numpy as np

'''
Sample script to retrieve data from NOAA CO-OPS FVCOM netcdf "file" (can be
OPeNDAP url), generate necessary grid topology (boundary info), and write 
GNOME compatible output.

'''
# specify local file or opendap url
data_dir = 'C:\\Users\\amy.macfadyen\\Documents\\Projects\\goods\\trunk\\static\\ocean_models\\PNNL'
#file1 = 'psm_0090.nc'
#filename = os.path.join(data_dir,file1)
files = ['psm_0090.nc','psm_0091.nc','psm_0092.nc','psm_0093.nc']
filenames = [os.path.join(data_dir,f) for f in files]

# the utools clasalish_sea requires a mapping of specific model variable names (values)
# to common names (keys) so that the clasalish_sea methods can work with FVCOM, SELFE,
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

# clasalish_sea instantiation creates a netCDF Dataset object as an attribute
salish_sea = tri_grid.ugrid(filenames)

# get longitude, latitude, and time variables
print 'Downloading data dimensions'
salish_sea.get_dimensions(var_map)
#fix time units
salish_sea.atts['time']['units'] = 'seconds since 2006-01-01 00:00:00'
#display available time range for model output
#nctools.show_tbounds(salish_sea.Dataset.variables['time'])
# UTM coordinates -- calculate lat/lon
x = salish_sea.data['lon'] #these are actually UTM even though named lon/lat in file
y = salish_sea.data['lat']
lon = np.ones_like(x); lat = np.ones_like(x)
for ii in range(len(x)):
    lat[ii], lon[ii] = nctools.utmToLatLng(10,x[ii],y[ii])
salish_sea.data['lon'] = lon
salish_sea.data['lat'] = lat

# get grid topo variables (nbe, nv)
print 'Downloading grid topo variables'
salish_sea.get_grid_topo(var_map)

# GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
salish_sea.atts['nbe']['order'] = 'cw'

# find and order the boundary
print 'Finding boundary'
bnd = salish_sea.find_bndry_segs()
print 'Ordering boundary'
seg_types = [0] * len(bnd)
salish_sea.order_boundary(bnd,seg_types)

# get the data
print 'Downloading data'
#salish_sea.get_data(var_map,tindex=[0,1,1]) #First time step only
salish_sea.get_data(var_map) #All time steps in file
 
print 'Writing to GNOME file'
salish_sea.write_unstruc_grid(os.path.join(data_dir, 'salish_sea_example.nc'))