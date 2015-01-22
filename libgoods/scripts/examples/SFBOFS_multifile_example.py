#!/usr/bin/env python
from libgoods import tri_grid, noaa_coops, data_files_dir, nctools
reload(tri_grid)
import datetime as dt
import os 
from netCDF4 import num2date, Dataset

'''
Sample script to retrieve data from unstructured grid netcdf "file" (can be
OPeNDAP url), generate necessary grid topology (boundary info), and write 
GNOME compatible output.

The boundary file is saved to the data files directory so it only needs 
to be generated once (unless you are subsetting the grid).

To process multiple files (urls) this script passes the filenames/urls in as a list 
-- this creates a netcdf4 MFDataset and is
a good option for not too many files (all output is written to one nc file for GNOME 
in this case)
Compare with sfbofs_multifile_example2 in this case, we use
a file list loop after the grid topo vars are loaded (as
this only has to be done once). 
'''
start = dt.date(2015,1,14)
end = dt.date(2015,1,21)
hour0 = 3
flist = noaa_coops.make_server_filelist('sfbofs',hour0,start,end=end,test_exist=False)

#generate text file with list of generated files
list_of_ofns = file(os.path.join(data_files_dir,'sfbofs_filelist.txt'), 'w')
list_of_ofns.write('NetCDF Files\n')

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

firsttime = 1

for f in flist[0:5]:
    
    if firsttime:
        
        firsttime = 0
        # class instantiation creates a netCDF Dataset object as an attribute -- 
        # use the first file in the list only
        sfbofs = tri_grid.ugrid(f)
        
#        # get longitude, latitude, and time variables
#        print 'Downloading data dimensions'
#        sfbofs.get_dimensions(var_map)
        
        # get grid topo variables (nbe, nv)
        print 'Downloading grid topo variables'
        sfbofs.get_grid_topo(var_map)
        # GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
        sfbofs.atts['nbe']['order'] = 'cw'
         
        # find and order the boundary
        print 'Finding boundary'
        bnd = sfbofs.find_bndry_segs()
        print 'Ordering boundary'
        seg_types = [0] * len(bnd)
        sfbofs.order_boundary(bnd,seg_types)
        
    else:
        
        sfbofs.update(f) 

    print 'Downloading data dimensions'
    sfbofs.get_dimensions(var_map)
    
    #get the data
    print 'Downloading data'
    #sfbofs.get_data(var_map,tindex=[0,1,1]) #First time step only
    sfbofs.get_data(var_map) #All time steps in file
    
    of_dt = nctools.round_time(num2date(sfbofs.data['time'][0],sfbofs.atts['time']['units']),roundto=3600)
    ofn = of_dt.strftime('%Y%m%d_%H') + '.nc'
    list_of_ofns.write('[FILE]  ' + ofn + '\n')
    print 'Writing to GNOME file'
    sfbofs.write_unstruc_grid(os.path.join(data_files_dir,ofn))
    
list_of_ofns.close()