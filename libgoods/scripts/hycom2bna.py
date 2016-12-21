#!/usr/bin/env python

"""
Script to read a HYCOM file and make a BNA map out of it for GNOME.

Usage:
  hycom2bna input_netcdf_file [outputfile]

(there are defaults in the script, but you probably want to pass in at
    least the input file name)

It draws a small land recangle for each grid box in the grid that is
does not have the water mask set.

This script was tested with one particular case of GNOME-compatilbe HYCONM files,
it may not work on others.
"""

from __future__ import print_function
import sys
import numpy as np
import netCDF4



try:
    infilename = sys.argv[1]
except IndexError:
    # use the default
    infilename = 'hycom_hindcast_currents_2006_surface_marianas_gnome.nc' 

try:
    outfilename = sys.argv[2]
except IndexError:
    # use the default
    outfilename = 'hycom_grid.bna' 

print("reading:", infilename)
nc = netCDF4.Dataset(infilename)


lat_var = nc.variables['lat']
lon_var = nc.variables['lon']
nx, ny = lat_var.shape

## add the extra row and column
lat = np.zeros((nx+1, ny+1))
lon = np.zeros((nx+1, ny+1))
lat[:-1,:-1] = lat_var[:,:]
lon[:-1,:-1] = lon_var[:,:]

# right column
lon[:,-1] = 2*lon[:,-2] - lon[:,-3]
lat[:,-1] = lat[:,-2]

# top row
lat[-1,:] = 2*lat[-2,:] - lat[-3,:]
lon[-1,:] = lon[-2,:]

# land-water mask
mask = nc.variables['mask'][:]

print("grid size is:",  nx, ny)

#generate the bna
print("writing out:", outfilename)
bna = open(outfilename, 'w')

# loop through the grid boxes:
count = 1
for i in range(nx):# nx
    for j in range(ny):# ny
        if not mask[i,j]:
            bna.write('"%i","1",5\n'%count)
            bna.write("%f, %f\n"%(lon[i,j], lat[i,j]) )
            bna.write("%f, %f\n"%(lon[i+1,j], lat[i+1,j]) )
            bna.write("%f, %f\n"%(lon[i+1,j+1], lat[i+1,j+1]) )
            bna.write("%f, %f\n"%(lon[i,j+1], lat[i,j+1]) )
            bna.write("%f, %f\n"%(lon[i,j], lat[i,j]) )
            count += 1

bna.close()

print("done")







