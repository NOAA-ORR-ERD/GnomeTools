#!/usr/bin/env python
from libgoods import utools, nctools, data_files_dir
import os 
import numpy as np

'''
Sample script to retrieve data from unstructured grid netcdf "file" (can be
OPeNDAP url), generate necessary grid topology (boundary info), and write 
GNOME compatible output.

The boundary file is saved to the data files directory so it only needs 
to be generated once (unless you are subsetting the grid).

To process multiple files (urls) either
a) pass the filenames/urls in as a list -- this creates a netcdf4 MFDataset and is
a good option for not too many files (all output is written to one nc file for GNOME 
in this case)
b) add a file list loop -- in this case put it after the grid topo vars are loaded (as
this only has to be done once). See NGOFS_multifile_example.py

'''

# specify local file or opendap url 
data_file = os.path.join(data_files_dir,'GalvCoarse_2010_0000.nc')

# the utools class requires a mapping of specific model variable names (values)
# to common names (keys) so that the class methods can work with FVCOM, SELFE,
# and ADCIRC which have different variable names
# (This seemed easier than finding them by CF long_names etc)
#!!!!!!!!txsuntans output on server does not include eles_surrounding_ele info
#I have it saved as a netcdf file included in libgoods data_files directory
var_map = { 
            'u_velocity':'uc', \
            'v_velocity':'vc', \
            'nodes_surrounding_ele':'cells',\
            'eles_surrounding_ele':'nbe',\
          }  

# class instantiation creates a netCDF Dataset object as an attribute
txsuntans = utools.ugrid(data_file)

# get longitude, latitude, and time variables
print 'Downloading data dimensions'
# Normally this would be a call to txself.get_dimensions(var_map) but
# this file has no time info, and only UTM coordinates -- so we do it
# manually

#Enter actual model time here for particular file
time_len = len(txsuntans.Dataset.dimensions['time'])
t_units = 'hours since 2014-01-01 00:00:00'
txsuntans.data['time'] = [t for t in range(time_len)]
txsuntans.atts['time'] = {'units':t_units}

x = txsuntans.Dataset.variables['xp'][:]
y = txsuntans.Dataset.variables['yp'][:]
lon = np.ones_like(x); lat = np.ones_like(x)
for ii in range(len(x)):
    lat[ii], lon[ii] = nctools.utmToLatLng(14,x[ii],y[ii])
txsuntans.data['lon'] = lon
txsuntans.data['lat'] = lat
txsuntans.atts['lon'] = {'long_name': 'longitude'}
txsuntans.atts['lat'] = {'long_name': 'latitude'}

# get grid topo variables (nbe, nv)
print 'Downloading grid topo variables'
txsuntans.get_grid_topo(var_map)


# GNOME requires boundary info -- this file can be read form data_files directory
# if saved or generated
print 'Loading/generating boundary segments'
bndry_file = os.path.join(data_files_dir, 'txsuntans.bry')
try:
    txsuntans.read_bndry_file(bndry_file)
except IOError:
    txsuntans.write_bndry_file('txsuntans',bndry_file)
    txsuntans.read_bndry_file(bndry_file)
txsuntans.data['nbe'] = txsuntans.data['nbe']
txsuntans.data['nv'] = txsuntans.data['nv']
# GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
txsuntans.atts['nbe']['order'] = 'ccw'

'''
!!!!!!!!!!!!!All the stuff above here only has to be done once -- if you want to 
process multiple files, I'd put a loop here and just keep overwriting 
txsuntans.data['u'] and ['v'] and incrementing the output file name
Also need to change txsuntans.data['time'] appropriately
'''
# get the data
print 'Loading u/v'
txsuntans.get_data(var_map,zindex=0) #First time step only

  
print 'Writing to GNOME file'
txsuntans.write_unstruc_grid(os.path.join(data_files_dir, 'txsuntans_example.nc'))