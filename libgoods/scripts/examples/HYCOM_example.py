from libgoods import curv_grid, nctools, data_files_dir
reload(curv_grid)
import numpy as np
import datetime
import os
from netCDF4 import num2date

url = 'http://tds.hycom.org/thredds/dodsC/glb_current_analysis'
sdate = datetime.datetime(2015,1,1,0,0)
edate = datetime.datetime(2015,1,5,0,0)

var_map = { 'time':'MT',
            'lon': 'Longitude',
            'lat': 'Latitude',
            'u': 'u',
            'v': 'v',
            }  
hycom = curv_grid.cgrid(url)
hycom.get_dimensions(var_map)
ts = num2date(hycom.data['time'][:],hycom.atts['time']['units'])
tid = np.where(np.logical_and(ts>=sdate,ts<=edate))

#adjust longitude: HYCOM lon goes from 74 --> 1019 ??
lon = np.mod(hycom.data['lon'],360)
hycom.data['lon'] = (lon > 180).choose(lon,lon-360)

#HYCOM uses time based on 1900 -- GNOME can't handle pre 1970, adjust before writing to file
native_tunits = hycom.atts['time']['units']
#hycom.data['time'],hycom.atts['time']['units'] = nctools.adjust_time(hycom.data['time'],hycom.atts['time']['units'])

hycom.subset([20,120,65,260]) #south lat, west lon, north lat, east lon

for num in range(len(tid)):
    print num
    hycom.get_data(var_map,tindex=[tid_start+num,tid_start+num+1,1],yindex=hycom.y,xindex=hycom.x,zindex=0,is3d=False)     
    hycom.data['time_ss'],hycom.atts['time']['units'] = nctools.adjust_time(hycom.data['time_ss'],native_tunits)    
    print num2date(hycom.data['time_ss'],hycom.atts['time']['units'])
    hycom.write_nc(os.path.join(data_files_dir,'HYCOM_example_ ' + str(num).zfill(3) + '.nc'),is3d=False)