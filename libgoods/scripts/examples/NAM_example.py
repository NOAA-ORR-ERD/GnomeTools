from __future__ import print_function

from libgoods import curv_grid, nctools, data_files_dir
import os


out_dir = data_files_dir
mon = 1
bbox = [36,-123,37.5,-121]

firsttime = True

for day in range(26,31):

    fl = 'http://data.oceansmap.com/thredds/dodsC/EDS/NAM5K/NAM5K2021' + str(mon).zfill(2) +str(day).zfill(2) + '.nc'
    print(fl)

    var_map = {'time':'time','longitude':'longitude','latitude':'latitude','u_velocity':'u','v_velocity':'v'}

    if firsttime:
        nam = curv_grid.cgrid(fl)
        nam.get_dimensions(var_map,get_time=False)
        nam.data['lon'] = nam.data['lon'] - 360
        nam.subset(bbox)
        firsttime = False

    nam.update(fl)

    nam.get_dimensions(var_map,get_xy=False)
    nam.get_data(var_map,tindex=[0,8,1],xindex=nam.x,yindex=nam.y)
    nam.data['time'] = nam.data['time'][0:8]
    nam.atts['wind'] = True
    nam.write_nc(os.path.join(out_dir,'NAM' + str(mon) + str(day) + '.nc'))
    
nctools.make_filelist_for_GNOME(out_dir,'NAM*.nc')