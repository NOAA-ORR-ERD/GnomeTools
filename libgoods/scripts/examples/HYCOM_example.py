from libgoods import curv_grid, nctools, data_files_dir
import numpy as np
import datetime
import os
from netCDF4 import num2date

#!!!!!!!!!
#These are things you might want to change
url = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.2'
sdate = datetime.datetime(2016,9,1,0,0) #First day to download
edate = datetime.datetime(2016,9,2,0,0) #Last day to download
#bbox = [25,-93.5,31.5,-85.5] #Geographic domain [South Lat, West Lon, North Lat, East Lon]
bbox = [
    55.5,
    175.,
    70.5,
    -160.
]

out_dir = data_files_dir
#out_dir = os.path.join(data_files_dir,'Indian_Ocean') #Where to write files (default is libgoods/data_files )
#!!!!!!!!


var_map = { 'time':'MT',
            'longitude': 'Longitude',
            'latitude': 'Latitude',
            'u_velocity': 'u',
            'v_velocity': 'v',
            'z': 'Depth',
            } 
            
hycom = curv_grid.cgrid(url)
hycom.get_dimensions(var_map)
ts = num2date(hycom.data['time'][:],hycom.atts['time']['units'])
tid = np.where(np.logical_and(ts>=sdate,ts<=edate))[0]
print 'Number of time steps:', len(tid)
#adjust longitude: HYCOM lon goes from 74 --> 1019 ??
lon = np.mod(hycom.data['lon'],360)
hycom.data['lon'] = (lon > 180).choose(lon,lon-360)

#HYCOM uses time based on 1900 -- GNOME can't handle pre 1970, adjust before writing to file
native_tunits = hycom.atts['time']['units']

#Determine geographic subset indices
hycom.subset(bbox,dl=1) #south lat, west lon, north lat, east lon
print hycom.x, hycom.y


for num,ti in enumerate(tid):
    print 'getting data time:', ti
    hycom.get_data(var_map,tindex=[ti,ti+1,1],yindex=hycom.y,xindex=hycom.x)     
    hycom.data['time_ss'],hycom.atts['time']['units'] = nctools.adjust_time(hycom.data['time_ss'],native_tunits)    
    print num2date(hycom.data['time_ss'],hycom.atts['time']['units'])
    hycom.write_nc(os.path.join(out_dir,'HYCOM_example_ ' + str(num).zfill(3) + '.nc'))
 
nctools.make_filelist_for_GNOME(out_dir,'HYCOM_example_*.nc',outfilename='HYCOM_filelist.txt')   