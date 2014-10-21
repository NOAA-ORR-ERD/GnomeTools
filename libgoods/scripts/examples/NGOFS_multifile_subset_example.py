#!/usr/bin/env python
from libgoods import utools, nctools, data_files_dir
import datetime as dt
import os 
import numpy as np

'''
Sample script to retrieve data from unstructured grid netcdf "file" (can be
OPeNDAP url), generate necessary grid topology (boundary info), and write 
GNOME compatible output.

The boundary file is saved to the data files directory so it only needs 
to be generated once (unless you are subsetting the grid).

To process multiple files (urls) we use
a file list loop after the grid topo vars are loaded (as
this only has to be done once). 
'''

sdate = dt.datetime(2013,3,1,3,0)
edate = dt.datetime(2013,3,7,3,0)

flist = []
url_stem = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/NGOFS/MODELS/'
while sdate < edate:
    yr = str(sdate.year)
    mon = str(sdate.month).zfill(2)
    day = str(sdate.day).zfill(2)
    hr = str(sdate.hour).zfill(2)
    url_sub_dir = yr + mon
    fname = 'nos.ngofs.fields.nowcast.' + yr + mon + day + '.t' + hr + 'z.nc'
    flist.append(url_stem + url_sub_dir + '/' + fname)
    delta_t = dt.timedelta(0.25,0)
    sdate = sdate + delta_t


# the utools class requires a mapping of specific model variable names (values)
# to common names (keys) so that the class methods can work with FVCOM, SELFE,
# and ADCIRC which have different variable names
# (This seemed easier than finding them by CF long_names etc)
var_map = { 'longitude':'lon', \
            'latitude':'lat', \
            'time':'time', \
            'u_velocity':'u', \
            'v_velocity':'v', \
            'nodes_surrounding_ele':'nv',\
            'eles_surrounding_ele':'nbe',\
          }  

# class instantiation creates a netCDF Dataset object as an attribute -- 
# use the first file in the list only
ngofs = utools.ugrid(flist[0])


# get longitude, latitude, and time variables
print 'Downloading data dimensions'
ngofs.get_dimensions(var_map)

#display available time range for model output
nctools.show_tbounds(ngofs.Dataset.variables['time'])

# get grid topo variables (nbe, nv)
print 'Downloading grid topo variables'
ngofs.get_grid_topo(var_map)

# GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
ngofs.atts['nbe']['order'] = 'cw'

#subsetting
print 'Subsetting'
nl = 31; sl = 28.5
wl = -89.25; el = -85

# Find all nodes and complete elements in subset box, lat/lon subset variables
ngofs.find_nodes_eles_in_ss(nl,sl,wl,el)

 # GNOME requires boundary info -- this file can be read form data_files directory
# if saved or generated
print 'Loading/generating boundary segments'
ss_bndry_file = os.path.join(data_files_dir, 'ngofs_mobile_sos.bry')
try:
    ngofs.read_bndry_file(ss_bndry_file)
except IOError:
    ngofs.write_bndry_file('ngofs',ss_bndry_file)
    ngofs.read_bndry_file(ss_bndry_file)
    
# get the data
ngofs.get_data(var_map,nindex=ngofs.nodes_in_ss)
print 'Downloading data'
for file in flist[1:]:
    print file
    ngofs_next = utools.ugrid(file)
    # check time units are consistent among files
    t = ngofs_next.Dataset.variables[var_map['time']]
    if t.units != ngofs.atts['time']['units']:
        print 'different time units among files'
        break
    ngofs_next.data['time'] = t[:]
    ngofs_next.data['nbe'] = ngofs.data['nbe']
    ngofs_next.eles_in_ss = ngofs.eles_in_ss
    ngofs_next.nodes_in_ss = ngofs.nodes_in_ss
    ngofs_next.get_data(var_map,nindex=ngofs.nodes_in_ss)
    # add data to original ngofs object
    ngofs.data['time'] = np.concatenate((ngofs.data['time'],ngofs_next.data['time']))
    ngofs.data['u'] = np.concatenate((ngofs.data['u'],ngofs_next.data['u']))
    ngofs.data['v'] = np.concatenate((ngofs.data['v'],ngofs_next.data['v']))
    

print 'Writing to GNOME file'
ngofs.write_unstruc_grid(os.path.join(data_files_dir, 'ngofs_multifile_ss_example.nc'))