from __future__ import print_function
from libgoods import curv_grid, nctools, noaa_coops, data_files_dir
reload(curv_grid)
import os
import datetime as dt

start = dt.date(2018,6,1)
end = dt.date(2018,6,3)
hour0 = 3
flist = noaa_coops.make_server_filelist('gomofs',hour0,start,end=end,test_exist=False)

var_map = { 'time':'ocean_time',
           }  
subset = 0
fnum = 0

gomofs = curv_grid.roms(flist[0])
gomofs.get_dimensions(var_map)
gomofs.data['lon'] = gomofs.data['lon_psi']
gomofs.data['lat'] = gomofs.data['lat_psi']
if subset:
    gomofs.subset([slat,wlon,nlat,elon])
    gomofs.get_grid_info(yindex=gomofs.y,xindex=gomofs.x)
else:
    gomofs.get_grid_info()


tlen = len(gomofs.data['time'])

for f in flist:
    print(f)    
    if fnum> 0:
        gomofs.update(f)
    
    if subset:
        gomofs.get_data(var_map,yindex=gomofs.y,xindex=gomofs.x)
    else:
        gomofs.get_data(var_map)
    fnum = fnum + 1
    
    outdir = os.path.join(data_files_dir,'gomofs_example')
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    ofn = os.path.join(outdir,'gomofs' + str(fnum).zfill(3) + '.nc')
    
    #gomofs.reduce_latlon_mesh_for_GNOME()
    gomofs.write_nc(ofn)
    
nctools.make_filelist_for_GNOME(outdir,'*.nc')
