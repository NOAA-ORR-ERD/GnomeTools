from libgoods import curv_grid, nctools, data_files_dir, noaa_coops
import datetime as dt
import os

flist = noaa_coops.make_server_filelist('nyofs',5,dt.date(2018,12,10),end=None,test_exist=True)
url = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NYOFS/fmrc/Aggregated_7_day_NYOFS_Fields_Forecast_best.ncd'

var_map = nctools.get_var_map(url,var_list=['longitude','latitude','time','u_velocity','v_velocity'])

nyofs = curv_grid.cgrid(url)
nyofs.get_dimensions(var_map)

subset = 1
tlen = len(nyofs.data['time'])

if subset:
    nyofs.subset([40.566,-74.114,40.654,-74.012])
    nyofs.get_grid_info(var_map,yindex=nyofs.y,xindex=nyofs.x)
    nyofs.get_data(var_map,yindex=nyofs.y,xindex=nyofs.x)
    ofn = os.path.join(data_files_dir,'nyofs_example.nc')
    nyofs.write_nc(ofn)
else:
    nyofs.get_grid_info(var_map)
    nyofs.get_data(var_map,tindex=[tlen-5,tlen,1])
    ofn = os.path.join(data_files_dir,'nyofs_example.nc')
    nyofs.write_nc(ofn)
