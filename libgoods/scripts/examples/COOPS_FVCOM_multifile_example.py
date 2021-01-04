#!/usr/bin/env python
from __future__ import print_function
from libgoods import tri_grid, noaa_coops, nctools
import datetime as dt
import os 
from netCDF4 import num2date
'''
Sample script to retrieve data from NOAA CO-OPS FVCOM netcdf "file" (can be
OPeNDAP url), generate necessary grid topology (boundary info), and write 
GNOME compatible output.

To script illustrates how to access data from multiple files (urls) by looping
through a filelist. Also only extracts a small subset from unstructured grid file

Alternatively, the list of filenames/urls can be passed directly when instantiating
the ugrid object-- this creates a netcdf4 MFDataset and isa good option for not too 
many files (all output is written to one nc file for GNOME in this case)

Since multiple files are created, also create a text file that can be loaded
into GNOME pointing to the individual files
'''

out_dir = 'ngofs'

# start = dt.date(2018,10,2)
# end = dt.date(2018,10,4)
# date = start
# dates = []
# while date <= end:
    # dates.append(date)
    # date += datetime.timedelta(days=1)	
#flist = ['fvcom_maine' + str(d.year) + str(d.month).zfill(2) + str(d.day).zfill(2) for d in dates]

flist = noaa_coops.make_server_filelist('ngofs',dt.date(2020,10,21))

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


# class instantiation creates a netCDF Dataset object as an attribute -- 
# use the first file in the list only
ngofs = tri_grid.ugrid(flist[0])

#get longitude, latitude
print('Downloading data dimensions')
ngofs.get_dimensions(var_map,get_time=False)

# get grid topo variables (nbe, nv)
print('Downloading grid topo variables')
ngofs.get_grid_topo(var_map)

 # subset bounding box
nl = 30.7; sl = 28.6
wl = -89.643; el = -87.391
ngofs.find_nodes_eles_in_ss(nl,sl,wl,el)

# find and order the boundary
print('Finding boundary')
bnd = ngofs.find_bndry_segs(subset=True)
print('Ordering boundary')
#In this case entire subset boundary will be set to land -- see COOPS_FVCOM_subset_example
#for how to use entire domain boundary to correctly determine type of subset boundary
ngofs.order_boundary(bnd)
# GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
ngofs.atts['nbe']['order'] = 'cw'

try:
    os.mkdir(out_dir)
except:
    pass

for f in flist:
    print(f)
    
    ngofs.update(f) 
    print('Downloading data dimensions')
    ngofs.get_dimensions(var_map,get_xy=False)
    of_dt = nctools.round_time(num2date(ngofs.data['time'][0],ngofs.atts['time']['units']),roundto=3600)
    ofn = of_dt.strftime('%Y%m%d_%H') + '.nc'
    
    if not os.path.exists(ofn):

        #get the data
        print('Downloading data')
        #ngofs.get_data(var_map) #First time step only

        ngofs.get_data(var_map,nindex=ngofs.nodes_in_ss) #All time steps in file
        

        print('Writing to GNOME file')
        ngofs.write_unstruc_grid(os.path.join(out_dir,ofn))
        
    else:
        print(ofn + ' aready exists')
    
nctools.make_filelist_for_GNOME(out_dir,'*.nc')
