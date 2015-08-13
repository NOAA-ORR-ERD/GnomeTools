#!/usr/bin/env python
from libgoods import curv_grid, nctools, noaa_coops, data_files_dir
import datetime as dt
import os 

# specify local file or opendap url

#sdate = dt.date(2015,3,8)
#edate = dt.date(2014,5,18)
#flist = noaa_coops.make_server_filelist('tbofs',0,sdate,test_exist=False)
#print "Number of files to aggregate: ", len(flist)

one_file = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/TBOFS/MODELS/201503/nos.tbofs.fields.n000.20150308.t00z.nc'


#We pick which variable we want to map to as this is sometimes not clear in virtual aggregations
var_map = { 'time':'ocean_time',
           }  
           
tbofs = curv_grid.roms(one_file)

tbofs.get_dimensions(var_map)
#nctools.show_tbounds(tbofs.Dataset.variables['ocean_time'])

# subset
nl = 28.2; sl = 27.67
wl = -82.5; el = -82.2
tbofs.subset([sl,wl,nl,el],lat='lat_psi',lon='lon_psi')

#tbofs.data['lon'] = tbofs.data['lon_psi']
#tbofs.data['lat'] = tbofs.data['lat_psi']

#get grid info
tbofs.get_grid_info(xindex=tbofs.x,yindex=tbofs.y)

print 'Getting data'
#get_data interps to rho grid unless tinterp=False
#tbofs.get_data(xindex=tbofs.x,yindex=tbofs.y)
tbofs.get_data(var_map,xindex=tbofs.x,yindex=tbofs.y) 

tbofs.data['lonc'] = tbofs.data['lon_rho']
tbofs.data['latc'] = tbofs.data['lat_rho']

    
print 'writing data'
#tbofs.write_nc(os.path.join(data_files_dir,'tbofs_example.nc'),is3d=False)