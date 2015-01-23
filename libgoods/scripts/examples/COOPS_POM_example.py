from libgoods import curv_grid, nctools, data_files_dir
import os


url = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/GBOFS/fmrc/Aggregated_7_day_GBOFS_Fields_Forecast_best.ncd'

var_map = nctools.get_var_map(url,var_list=['lon','lat','time','u','v'])

gbofs = curv_grid.cgrid(url)
gbofs.get_dimensions(var_map)

subset = 1
tlen = len(gbofs.data['time'])

if subset:
    gbofs.subset([29.104,-95.124,29.454,-94.613])
    gbofs.get_grid_info(yindex=gbofs.y,xindex=gbofs.x)
    gbofs.get_data(var_map,yindex=gbofs.y,xindex=gbofs.x)
    ofn = os.path.join(data_files_dir,'gbofs_example.nc')
    gbofs.write_nc(ofn)
else:
    gbofs.get_grid_info()
    gbofs.get_data(var_map,tindex=[tlen-5,tlen,1])
    ofn = os.path.join(data_files_dir,'gbofs_example.nc')
    gbofs.write_nc(ofn)
