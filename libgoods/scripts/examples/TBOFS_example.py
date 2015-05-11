#!/usr/bin/env python
from libgoods import curv_grid, nctools, noaa_coops, data_files_dir
import datetime as dt
import os 

# specify local file or opendap url
#data_url = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NGOFS/MODELS/201309/nos.ngofs.fields.f000.20130901.t03z.nc'
sdate = dt.date(2015,3,8)
#edate = dt.date(2014,5,18)

flist = noaa_coops.make_server_filelist('tbofs',0,sdate,test_exist=False)
print "Number of files to aggregate: ", len(flist)

#We pick which variable we want to map to as this is sometimes not clear in virtual aggregations
var_map = { 'time':'ocean_time',
           }  
           
tbofs = curv_grid.roms(flist) #just do the first couple for testing
tbofs.get_dimensions(var_map)
nctools.show_tbounds(tbofs.Dataset.variables['ocean_time'])

#get grid info
tbofs.get_grid_info()

# subset
#nl = 28.2; sl = 27.67
#wl = -82.9; el = -82.2
#tbofs.subset(nl,sl,wl,el)

print 'Getting data'
#get_data interps to rho grid unless tinterp=False
#tbofs.get_data(xindex=tbofs.x,yindex=tbofs.y)
tbofs.get_data(var_map) 

# GNOME needs lat/lon mesh reduced by last row/column -- lat/lon on psi grid, u/v on rho grid
tbofs.reduce_latlon_mesh_for_GNOME()
ofn = os.path.join(data_files_dir,'tbofs_rho_reduced.nc')
tbofs.data['lon'] = tbofs.data['lon_psi']
tbofs.data['lat'] = tbofs.data['lat_psi']

    
print 'writing data'
tbofs.write_nc(os.path.join(data_files_dir,'tbofs_example.nc'),is3d=False)