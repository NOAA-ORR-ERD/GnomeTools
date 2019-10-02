#!/usr/bin/env python
from __future__ import print_function
from libgoods import tri_grid, noaa_coops, data_files_dir, nctools
import datetime as dt
import os 
from netCDF4 import num2date

'''
Sample script to retrieve data from NOAA CO-OPS FVCOM netcdf "file" (can be
OPeNDAP url), generate necessary grid topology (boundary info), and write 
GNOME compatible output.

To script illustrates how to access data from multiple files (urls) by looping
through a filelist.

Alternatively, the list of filenames/urls can be passed directly when instantiating
the ugrid object-- this creates a netcdf4 MFDataset and isa good option for not too 
many files (all output is written to one nc file for GNOME in this case)

Since multiple files are created, also create a text file that can be loaded
into GNOME pointing to the individual files
'''
start = dt.date(2019,3,1)
end = dt.date(2019,3,8)
hour0 = 3
flist = noaa_coops.make_server_filelist('creofs',hour0,start,end=end,test_exist=False)

#generate text file with list of generated files
list_of_ofns = file(os.path.join(data_files_dir,'creofs_filelist.txt'), 'w')
list_of_ofns.write('NetCDF Files\n')

nl = 47; sl = 45
wl = -126.2; el = -123.6

# the utools class requires a mapping of specific model variable names (values)
# to common names (keys) so that the class methods can work with FVCOM, SELFE,
# and ADCIRC which have different variable names
# (This seemed easier than finding them by CF long_names etc)
var_map = { 'longitude':'lon', \
            'latitude':'lat', \
            'time':'time', \
            'u_velocity':'u', \
            'v_velocity':'v', \
            'nodes_surrounding_ele':'ele',\
            'eles_surrounding_ele':'',\
          }  

firsttime = 1

for f in flist:
    
    if firsttime:
        
        firsttime = 0
        # class instantiation creates a netCDF Dataset object as an attribute -- 
        # use the first file in the list only
        creofs = tri_grid.ugrid(f)
        
#        # get longitude, latitude, and time variables
#        print 'Downloading data dimensions'
#        creofs.get_dimensions(var_map)
        
        # get grid topo variables (nbe, nv)
        print('Downloading grid topo variables')
        creofs.get_grid_topo(var_map)
        creofs.get_dimensions(var_map,get_time=False)
        # GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
        creofs.atts['nbe']['order'] = 'ccw'
        
        creofs.find_nodes_eles_in_ss(nl,sl,wl,el)
        
        # find and order the boundary
        print('Finding boundary')
        bnd = creofs.find_bndry_segs(subset=True)
        #In order to correctly specify land/ow segments requires comparison with full domain boundary
        #Create this by downloading entire domain grid info then saving it (write_bndry_file)
        bry_file = 'C:\\Users\\amy.macfadyen\\PyProjects\\goods\\trunk\\static\\ocean_models\\COOPS\\creofs.bry'
        land_nodes = creofs.find_subset_land_nodes(bry_file)
        #seg_types = noaa_coops.specify_bnd_types('creofs',bnd,ss_land_nodes=land_nodes)
        print('Ordering boundary')
        creofs.order_boundary(bnd)
        
        out_dir = os.path.join(data_files_dir,'creofs')
        try:
            os.mkdir(out_dir)
        except:
            pass
        
    else:
        
        creofs.update(f) 

    print('Downloading data dimensions')
    creofs.get_dimensions(var_map,get_xy=False)
    
    #get the data
    print('Downloading data')
    #creofs.get_data(var_map,tindex=[0,1,1]) #First time step only
    creofs.get_data(var_map,zindex=-1,nindex=creofs.nodes_in_ss) #All time steps in file
    
    of_dt = nctools.round_time(num2date(creofs.data['time'][0],creofs.atts['time']['units']),roundto=3600)
    ofn = of_dt.strftime('%Y%m%d_%H') + '.nc'
    list_of_ofns.write('[FILE]  ' + ofn + '\n')
    print('Writing to GNOME file')
    creofs.write_unstruc_grid(os.path.join(out_dir,ofn))
    
list_of_ofns.close()
