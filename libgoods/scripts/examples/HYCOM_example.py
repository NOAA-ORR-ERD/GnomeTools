from __future__ import print_function
from libgoods import rect_model, nctools, data_files_dir
import numpy as np
import cftime
import os
from netCDF4 import num2date
#reload(rect_model)
#!!!!!!!!!
#These are things you might want to change
#url = 'http://eds.ioos.us/thredds/dodsC/ioos/hycom_navy/Navy_HYCOM_best.ncd'
url = 'https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0'
#url = 'http://tds.hycom.org/thredds/dodsC/GLBy0.08/latest'
sdate = cftime.num2date(0,'hours since 2020-07-10 00:00:00.000 UTC') #First day to download
edate = cftime.num2date(0,'hours since 2020-07-20 00:00:00.000 UTC') #Last day to download
#bbox = [25,-93.5,31.5,-85.5] #Geographic domain [South Lat, West Lon, North Lat, East Lon]
bbox = [
    55,
    170.0,
    67.5,
    -160.0
]

out_dir = data_files_dir



var_map = { 'time':'time',
            'lon': 'lon',
            'lat': 'lat',
            'u_velocity': 'water_u',
            'v_velocity': 'water_v',
            'z': 'Depth',
            } 
            
hycom = rect_model.rect(url)
hycom.get_dimensions(var_map)
ts = cftime.num2date(hycom.time[:],hycom.time_units)
tid = np.where(np.logical_and(ts>=sdate,ts<=edate))[0]
print('Number of time steps:', len(tid))
print([tid[0],tid[-1],1])

hycom.subset(bbox,dl=1) #south lat, west lon, north lat, east lon
print(hycom.x, hycom.y)

print('writing file')
ofn = os.path.join(out_dir,'HYCOM_0710-0720.nc')
hycom.write_nc(var_map,ofn,t_index=[tid[0],tid[-1],1],dl=1)

 

