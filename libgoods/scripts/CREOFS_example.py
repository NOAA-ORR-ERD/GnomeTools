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

To process multiple files (urls) either
a) pass the filenames/urls in as a list -- this creates a netcdf4 MFDataset and is
a good option for not too many files (all output is written to one nc file for GNOME 
in this case)
b) add a file list loop -- in this case put it after the grid topo vars are loaded (as
this only has to be done once). See NGOFS_multifile_example.py

'''


# specify local file or opendap url -- in this case files are one time step, not aggregated
data_url = 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/CREOFS/MODELS/201308/nos.creofs.fields.f000.20130820.t03z.nc'
'''Some CO-OPS model notes:
nos.OFS.fields.fHHH.YYYYMMDD.tCCz.nc
OFS  is cbofs, dbofs, tbofs, ngofs, or creofs to stand for the name of operational forecast systems for Chesapeake Bay, Delaware Bay, Tampa Bay, North Gulf of Mexico, and Columbia River Estuary.
CC is the model cycle runtime (i.e. 00, 06, 12, 18 for CBOFS, DBOFS, and TBOFS, 03,09,15,21 for NGOFS and CREOFS)
HH is 2 digits for the nowcast|forecast hour of the product from 00 - 48
HHH is 3 digits for the nowcast|forecast hour of the product from 000 - 048
YYYY is the calendar year of model runtime
MM is the calendar month of model runtime
DD is the calendar date of model runtime'''

# generate a file list of all files for this forecast
fstem = data_url.split('nos.creofs.fields')[0]
fname = data_url.split('/')[-1].split('f000')
flist = []
for hh in range(49):
    flist.append(fstem + fname[0] + 'f' + str(hh).zfill(3) + fname[1])

# the utools class requires a mapping of specific model variable names (values)
# to common names (keys) so that the class methods can work with FVCOM, SELFE,
# and ADCIRC which have different variable names
# (This seemed easier than finding them by CF long_names etc)
#!!!!!!!!CREOFS output on server does not include eles_surrounding_ele info
#I have it saved as a netcdf file included in libgoods data_files directory
var_map = { 'longitude':'lon', \
            'latitude':'lat', \
            'time':'time', \
            'u_velocity':'u', \
            'v_velocity':'v', \
            'nodes_surrounding_ele':'ele',\
            'eles_surrounding_ele':'',\
          }  

# class instantiation creates a netCDF Dataset object as an attribute
# creofs = utools.ugrid(flist) #multiple files
creofs = utools.ugrid(data_url) #single file output

# get longitude, latitude, and time variables
print 'Downloading data dimensions'
creofs.get_dimensions(var_map)

#display available time range for model output
nctools.show_tbounds(creofs.time)

# get grid topo variables (nbe, nv)
print 'Downloading grid topo variables'
try:
    creofs.get_grid_topo(var_map)
except KeyError: #model output on server doesn't have nbe
    creofs.build_face_face_connectivity()
# GNOME needs to know whether the elements are ordered clockwise (FVCOM) or counter-clockwise (SELFE)
creofs.atts['nbe']['order'] = 'ccw'

# GNOME requires boundary info -- this file can be read form data_files directory
# if saved or generated
print 'Loading/generating boundary segments'
bndry_file = os.path.join(data_files_dir, 'creofs.bry')
try:
    creofs.read_bndry_file(bndry_file)
except IOError:
    creofs.write_bndry_file('creofs',bndry_file)
    creofs.read_bndry_file(bndry_file)

# get the data
print 'Downloading data'
#creofs.get_data(var_map,tindex=[0,1,1]) #First time step only
creofs.get_data(var_map,zindex=-1) #All time steps in file
 
print 'Writing to GNOME file'
creofs.write_unstruc_grid(os.path.join(data_files_dir, 'creofs_example.nc'))