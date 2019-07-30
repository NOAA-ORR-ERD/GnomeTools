#!/usr/bin/env python
from __future__ import print_function
from libgoods import curv_grid, nctools, noaa_coops, data_files_dir
reload(curv_grid)
import datetime as dt
import os 

# specify local file or opendap url

# sdate = dt.date(2019,6,10)
# edate = dt.date(2019,6,11)
# flist = noaa_coops.make_server_filelist('tbofs',0,sdate,test_exist=False)
# print("Number of files to aggregate: ", len(flist))

# flist = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/TBOFS/MODELS/201906/nos.tbofs.fields.n006.20190605.t12z.nc'
flist = 'C:\\Users\\amy.macfadyen\\Documents\\COOPS_tests\\nos.tbofs.fields.n006.20190605.t12z.nc'

#We pick which variable we want to map to as this is sometimes not clear in virtual aggregations
var_map = { 'time':'ocean_time',
           }  
           
tbofs = curv_grid.roms(flist)

print('getting dimensions')
tbofs.get_dimensions(var_map)
#nctools.show_tbounds(tbofs.Dataset.variables['ocean_time'])

print('getting grid info')
#get grid info
tbofs.get_grid_info(is3d=False)

print('Getting data')
#get_data interps to rho grid unless interp=False
tbofs.get_data(var_map,interp=False) 

print('writing data')
tbofs.write_nc_native(os.path.join(data_files_dir,'tbofs.nc'))
