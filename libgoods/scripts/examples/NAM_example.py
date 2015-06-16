#!/usr/bin/env python
from libgoods import curv_grid, nctools, data_files_dir
import datetime as dt
import os 
import pyproj

url = 'http://thredds.ucar.edu/thredds/dodsC/grib/NCEP/NAM/CONUS_12km/Best'
var_map = { 'time':'time1',
            'lon': 'x',
            'lat': 'y',
            'u': 'u-component_of_wind_height_above_ground',
            'v': 'v-component_of_wind_height_above_ground',
            }  
            
nam = curv_grid.cgrid(url)
nam.get_dimensions(var_map)

#convert x,y in lambert conformal conic to lat/lon
p = Proj(r'+proj=lcc +lat_0=25.0 +lon_0=265) #N


latitude_of_projection_origin: 25.0
longitude_of_central_meridian: 265.0
standard_parallel: 25.0
earth_radius: 6371229.0