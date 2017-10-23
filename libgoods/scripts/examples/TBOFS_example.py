#!/usr/bin/env python
from __future__ import print_function
from libgoods import curv_grid, nctools, noaa_coops, data_files_dir
import datetime as dt
import os 

# specify local file or opendap url

#sdate = dt.date(2015,3,8)
#edate = dt.date(2014,5,18)
#flist = noaa_coops.make_server_filelist('tbofs',0,sdate,test_exist=False)
#print "Number of files to aggregate: ", len(flist)

one_file = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/TBOFS/MODELS/201703/nos.tbofs.fields.n000.20170308.t00z.nc'


#We pick which variable we want to map to as this is sometimes not clear in virtual aggregations
var_map = { 'time':'ocean_time',
            'z':'s_rho',
           }  
           
tbofs = curv_grid.roms(one_file)

tbofs.get_dimensions(var_map)
#nctools.show_tbounds(tbofs.Dataset.variables['ocean_time'])


#get grid info
tbofs.get_grid_info(is3d=True)

print('Getting data')
#get_data interps to rho grid unless tinterp=False

tbofs.get_data(var_map) 

print('writing data')
tbofs.write_nc(os.path.join(data_files_dir,'tbofs_3d.nc'),is3d=True)
