#!/usr/bin/env python

# code to write GNOME compatible netcdf files.

import numpy as np

from netCDF4 import Dataset

def write_curv_grid(time, lon, lat, u, v, atts, ofn):
    """

    Write GNOME compatible netCDF file (netCDF3) from regular grid data
     Call signature::
       subset(time,lon,lat,u,v,atts,dfn)
 
    Arguments:
      *time*: 1-D time vector
      *lon* : 2-D in decimal degrees -180:180
      *lat* : 2-D in decimal degrees -90:90
      *u*   : 3-D eastward velocity component
      *v*   : 3-D northward velocity component
      *atts*: nested dict object of variable attributes with keys u, v, t, [mask]
      *ofn* : netCDF output file name -- string

    """
    nc = Dataset(ofn,'w',format='NETCDF3_CLASSIC')
    # Global Attributes
    #setattr(nc,'grid_type','curvilinear')
    nc.grid_type = 'curvilinear'

    # add Dimensions
    nc.createDimension('x',lon.shape[1])
    nc.createDimension('y',lat.shape[0])
    nc.createDimension('time', None)

    # create variables
    nc_time = nc.createVariable('time','f4',('time',))
    nc_lon = nc.createVariable('lon','f4',('y','x'))
    nc_lat = nc.createVariable('lat','f4',('y','x'))
    nc_u = nc.createVariable('water_u','f4',('time','y','x'))
    nc_v = nc.createVariable('water_v','f4',('time','y','x'))
    # add data
    nc_lon[:] = lon
    nc_lat[:] = lat
    if len(u.shape) == 2:
        nc_time[0] = time
        nc_u[0,:] = u
        nc_v[0,:] = v
    else:
        nc_time[:] = time
        nc_u[:] = u
        nc_v[:] = v
    if atts.has_key('mask'):
        nc_mask = nc.createVariable('mask','f4',('y','x'))
        nc_mask[:] = atts['mask']
    # add variable attributes from 'atts' (nested dict object)
    for an_att in atts['t'].iteritems():
        setattr(nc_time,an_att[0],an_att[1])

    for att, val in atts['u'].iteritems():
        setattr(nc_u, att, val)

    for an_att in atts['v'].iteritems():
        setattr(nc_v,an_att[0],an_att[1])

    nc.close()


if __name__=="__main__":
    # create some fake data
    time = np.arange(10).astype(np.float) / 24.0
    nt = len(time)
    print time
    nx, ny = 4, 6
    
    lon = np.random.random((ny,nx)) * -30 - 40 
    lat = np.random.random((ny,nx)) * 20 + 20 
    u = np.random.random((nt,ny,nx)) * 4 - 1
    v = np.random.random((nt,ny,nx)) * 4 - 1
    
    u_atts = {'long_name': "Eastward Water Velocity",
              'units':"m/s",
              'missing_value':"-99.f",
              'FillValue':   "-99.f",
              'standard_name': "eastward_sea_water_velocity",
              }
    v_atts = {'long_name': "Northward Water Velocity",
              'units':"m/s",
              'missing_value':"-99.f",
              'FillValue':   "-99.f",
              'standard_name': "northward_sea_water_velocity",
              }
    t_atts = {'long_name':"Time",
              'base_date': "1995, 1, 1, 0",
              'units': "days since 1995-01-01 0:00:00 00:00",
              'standard_name': "time",
              }
    atts = {'u': u_atts,
            'v': v_atts,
            't': t_atts,
            }
    
    write_curv_grid(time, lon, lat, u, v, atts, 'test.nc')


