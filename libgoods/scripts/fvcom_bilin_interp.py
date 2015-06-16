# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 16:07:23 2015

@author: amy.macfadyen
"""
from libgoods import data_files_dir
from netCDF4 import Dataset
import os
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

nc = Dataset(os.path.join(data_files_dir,'ngofs_ss_example.nc'))
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
lonc = nc.variables['lonc'][:]
latc = nc.variables['latc'][:]
nv = nc.variables['nv'][:]
nbe = nc.variables['nbe'][:]
u = nc.variables['u'][:]
v = nc.variables['v'][:]

tri = mtri.Triangulation(lon, lat, (nv-1).transpose())
plt.triplot(tri, 'k.-')

plt.draw()
plt.show()