#!/usr/bin/env python
from libgoods import tri_grid, nctools, noaa_coops, data_files_dir
import os 

'''
Sample script to retrieve data from NOAA CO-OPS FVCOM netcdf "file" (can be
OPeNDAP url), generate necessary grid topology (boundary info), and write 
GNOME compatible output.

'''
# specify local file or opendap url
data_url = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NGOFS/MODELS/201501/nos.ngofs.fields.n000.20150115.t09z.nc'
#data_url = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/ngofs/MODELS/201501/nos.ngofs.fields.f001.20150102.t15z.nc'
# the utools class requires a mapping of specific model variable names (values)
# to common names (keys) so that the class methods can work with FVCOM, SELFE,
# and ADCIRC which have different variable names
# (This seemed easier than finding them by CF long_names etc)
var_map = { 'longitude':'lon', \
            'latitude':'lat', \
            'lonc': 'lonc',\
            'latc': 'latc',\
            'a1u': 'a1u',\
            'a2u': 'a2u',\
            'time':'time', \
            'u_velocity':'u', \
            'v_velocity':'v', \
            'nodes_surrounding_ele':'nv',\
            'eles_surrounding_ele':'nbe',\
          }  

# class instantiation creates a netCDF Dataset object as an attribute
ngofs = tri_grid.ugrid(data_url)

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

nl = 29.3; sl = 28.7
wl = -90; el = -89
ngofs.find_nodes_eles_in_ss(nl,sl,wl,el)

# find and order the boundary
print 'Finding boundary'
bnd = ngofs.find_bndry_segs(subset=True)
#In order to correctly specify land/ow segments requires comparison with full domain boundary
#Create this by downloading entire domain grid info then saving it (write_bndry_file)
bry_file = 'C:\\Users\\amy.macfadyen\\Documents\\Projects\\goods\\trunk\\static\\ocean_models\\COOPS\\ngofs.bry'
land_nodes = ngofs.find_subset_land_nodes(bry_file)
seg_types = noaa_coops.specify_bnd_types('ngofs',bnd,ss_land_nodes=land_nodes)
print 'Ordering boundary'
ngofs.order_boundary(bnd,seg_types)

# get the data
print 'Downloading data'
#ngofs.get_data(var_map,tindex=[0,1,1]) #First time step only
ngofs.get_data(var_map,nindex=ngofs.nodes_in_ss) #All time steps in file
 
print 'Writing to GNOME file'
ngofs.write_unstruc_grid(os.path.join(data_files_dir, 'ngofs_ss_example.nc'))