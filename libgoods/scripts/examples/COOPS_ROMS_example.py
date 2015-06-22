from libgoods import curv_grid, nctools, noaa_coops, data_files_dir
import os
import datetime as dt


start = dt.date(2015,3,9)
end = dt.date(2015,3,9)
hour0 = 0
flist = noaa_coops.make_server_filelist('dbofs',hour0,start,end=end,test_exist=False)

var_map = { 'time':'ocean_time',
           }  
subset = 0
slat = 38.8
nlat = 39.2
wlon = -75.4
elon = -74.8

fnum = 0

dbofs = curv_grid.roms(flist[0])
dbofs.get_dimensions(var_map)
dbofs.data['lon'] = dbofs.data['lon_psi']
dbofs.data['lat'] = dbofs.data['lat_psi']
if subset:
    dbofs.subset([slat,wlon,nlat,elon])
    dbofs.get_grid_info(yindex=dbofs.y,xindex=dbofs.x)
else:
    dbofs.get_grid_info()


tlen = len(dbofs.data['time'])



for f in flist:
    print f    
    if fnum> 0:
        dbofs.update(f)
    
    if subset:
        dbofs.get_data(var_map,yindex=dbofs.y,xindex=dbofs.x)
    else:
        dbofs.get_data(var_map)
    fnum = fnum + 1
    
    outdir = os.path.join(data_files_dir,'dbofs_example')
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    ofn = os.path.join(outdir,'dbofs' + str(fnum).zfill(3) + '.nc')
    
    dbofs.reduce_latlon_mesh_for_GNOME()
    dbofs.write_nc(ofn)
    
nctools.make_filelist_for_GNOME(outdir,'*.nc')