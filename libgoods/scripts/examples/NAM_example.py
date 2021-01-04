from __future__ import print_function
from libgoods import curv_model, nctools, data_files_dir
import os
import datetime
'''
This accesses NAM from old IOOS EDS server (isn't on new server yet 08/2020)
You may need to change the stem in the filename below based on which NAM
is desired -- NAMAK, NAMHI, NAM5K etc.
See http://data.oceansmap.com/eds_thredds/catalog/EDS/catalog.html to
determine what's available
'''
sdate = datetime.date(2020,7,20)
edate = datetime.date(2020,8,2)
#edate = datetime.date(2020,8,2)

out_dir = data_files_dir
bbox = [
    61.69,
    177.31,
    65.9,
    -163.353
]

var_map = {'time':'time','lon':'longitude','lat':'latitude','u_velocity':'u','v_velocity':'v'}

firsttime = True
numdays = (edate - sdate).days

date_list = [sdate + datetime.timedelta(days = t) for t in range(0,numdays)]

for this_date in date_list:

    fl = 'http://data.oceansmap.com/eds_thredds/dodsC/EDS/NAMAK/NAMAK2020' + str(this_date.month).zfill(2) +str(this_date.day).zfill(2) + '.nc'
    print(fl)

    if firsttime:
        nam = curv_model.curv(fl)
        nam.get_dimensions(var_map,get_time=False)
        nam.subset(bbox,dl=1)
        firsttime = False

    nam.update(fl)

    nam.get_dimensions(var_map,get_xy=False)
    ofn = os.path.join(out_dir,'NAM_' + this_date.strftime("%Y_%m%d") +'.nc')
    nam.write_nc_gnome1(var_map,ofn,t_index=[0,8,1],dl=True,wind=True)

    
nctools.make_filelist_for_GNOME(out_dir,'NAM*.nc')