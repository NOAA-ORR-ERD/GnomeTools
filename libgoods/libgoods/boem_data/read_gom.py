#!/usr/bin/env python

"""
Code to read times from MMS/Oey GOM binary format

It also includes code to write it out in a GNOME-compatible netcdf format.

"""

import datetime
import numpy as np

print "importing netCDF4"
from netCDF4 import Dataset
print "done importing netCDF4"


class data_reader():

    time_dt = np.dtype([('year', '>i4'),('month','>i4'),('day','>i4'),('hour','>i4')])
    data_dt = np.dtype('>f4')
    #np.dtype([('uo', '>f4'),('vo','>f4'),('uw','>f4'),('vw','>f4'),('el','>f4')])
    fields = {'wind_u':    0,
              'wind_v':    1,
              'ocean_u':   2,
              'ocean_v':   3,
              'elevation': 4,
              }
    def __init__(self, infilename, grid_size):
        self.infilename = infilename

        self.grid_size = grid_size
        if grid_size == 'coarse':
            self.nx, self.ny = (200, 200)
            self.num_timesteps = 480
        elif grid_size == 'fine':
            self.nx, self.ny = (320, 428)
            self.num_timesteps = 240
        else:
            raise ValueError('grid_size must be "coarse" or "fine"')
    
    def get_field(self, field_name):
        i = self.fields[field_name]
        return self.data[:,i,:,:]
    
    def read_data(self):
        print "reading data from:", self.infilename
        infile = file(self.infilename, 'rb')

        self.data = np.empty((self.num_timesteps, 5, self.nx, self.ny), dtype=self.data_dt)

        # read the data:
        for ts in range(self.num_timesteps):
            if ts%20 == 0: print "reading time step:", ts
            #read the time stamp
            time = np.fromfile(infile, count=1, dtype=self.time_dt)
            # load the data:
            for i in range(5):
                self.data[ts, i, :] = np.fromfile( infile, count=(self.nx*self.ny), dtype=self.data_dt ).reshape((self.nx, self.ny), order='F')
        infile.close()
        print "done reading data block"

    def read_times(self):

        infile = file(self.infilename, 'rb')

        times = []
        for ts in range(self.num_timesteps):
            #read the time stamp
            time = np.fromfile(infile, count=1, dtype=self.time_dt)
            time['year'] += 1900
            if time['hour'] == 24:
                # some of  the time stamps have hour=24
                # I'm assuming that that should really be the 0th hour of the next day.
                # this code should fix that.
                time['hour'] -= 1 # to get around the "24th hour" problem
                dt = datetime.datetime(*time[0])
                dt = dt + datetime.timedelta(hours=1)
            else:
                dt = datetime.datetime(*time[0])
            times.append(dt)

            ##skip the data
            infile.seek(self.data_dt.itemsize*5*self.nx*self.ny, 1)
        infile.close()
        self.times = times

    def read_grid(self):
        infile = file(self.infilename, 'rb')
        
        #Amy comment: this is created but never used?
        #self.grid = np.empty((self.nx, self.ny, 3), dtype=self.data_dt)

        # skip the data
        bytes_offset = self.num_timesteps * (self.time_dt.itemsize + self.data_dt.itemsize * 5 * self.nx * self.ny)
        infile.seek(bytes_offset, 0)

        # read the data:
        infile.read(16)
        # read the grid
        self.lon = np.fromfile(infile,
                           count = self.nx*self.ny,
                           dtype=self.data_dt).reshape((self.nx, self.ny), order='F')
        self.lat = np.fromfile(infile,
                           count = self.nx*self.ny,
                           dtype=self.data_dt).reshape((self.nx, self.ny), order='F')
        self.depth = np.fromfile(infile,
                           count = self.nx*self.ny,
                           dtype=self.data_dt).reshape((self.nx, self.ny), order='F')
        infile.close()
    
    def read_all(self):
        self.read_times()
        self.read_data()
        self.read_grid()

    def to_netcdf(self, filename):
        """
        takes a MMS data class, and converts to format netcdf wants
        """
        ref_time = datetime.datetime(1993, 1, 1, 0)
        print "ref_time", ref_time
        
        # convert time to days since start time
        time = np.array([TD2float_days(t-ref_time) for t in self.times], dtype=np.float32)
        
        # lat-lon (ny, nx) -- lat first!
        lat = self.lat.transpose()
        lon = self.lon.transpose()

        # velocities (Water)
        u = self.get_field('ocean_u').transpose((0,2,1))
        # velocities (Water)
        v = self.get_field('ocean_v').transpose((0,2,1))
        
        base_date = "%i, %i, %i, %i"%(ref_time.year, ref_time.month, ref_time.day, ref_time.hour)
        # fixme: this isn't doing hours!
        time_units = "days since %i-%02i-%02i 0:00:00 00:00"% (ref_time.year, ref_time.month, ref_time.day)
        
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
                  'base_date': base_date,
                  'units': time_units,
                  'standard_name': "time",
                  }
        atts = {'u': u_atts,
                'v': v_atts,
                't': t_atts,
                }
    
        print "writing netcdf file"

        self.write_curv_grid(time, lon, lat, u, v, atts, filename)
        
    def to_netcdf_wind(self, filename):
            """
            takes a MMS data class, and converts to format netcdf wants
            """
            ref_time = datetime.datetime(1993, 1, 1, 0)
            print "ref_time", ref_time
            
            # convert time to days since start time
            time = np.array([TD2float_days(t-ref_time) for t in self.times], dtype=np.float32)
            
            # lat-lon (ny, nx) -- lat first!
            lat = self.lat.transpose()
            lon = self.lon.transpose()
    
            # velocities (Water)
            u = self.get_field('wind_u').transpose((0,2,1))
            # velocities (Water)
            v = self.get_field('wind_v').transpose((0,2,1))
            
            base_date = "%i, %i, %i, %i"%(ref_time.year, ref_time.month, ref_time.day, ref_time.hour)
            # fixme: this isn't doing hours!
            time_units = "days since %i-%02i-%02i 0:00:00 00:00"% (ref_time.year, ref_time.month, ref_time.day)
            
            u_atts = {'long_name': "Eastward Surface Air Velocity",
                      'units':"m/s",
                      'missing_value':"-99.f",
                      'FillValue':   "-99.f",
                      'standard_name': "eastward_surface_air_velocity",
                      }
            v_atts = {'long_name': "Northward Surface Air Velocity",
                      'units':"m/s",
                      'missing_value':"-99.f",
                      'FillValue':   "-99.f",
                      'standard_name': "northward_surface_air_velocity",
                      }
            t_atts = {'long_name':"Time",
                      'base_date': base_date,
                      'units': time_units,
                      'standard_name': "time",
                      }
            atts = {'u': u_atts,
                    'v': v_atts,
                    't': t_atts,
                    }
        
            print "writing netcdf file"
    
            self.write_curv_grid_wind(time, lon, lat, u, v, atts, filename)
    

    def subset_data():
        """
        take a subset of the data
        """
        min_x, max_x = 75, 125
        min_y, max_y = 75, 125
        min_time, max_time = 100, 120

        self.times = self.times[min_time, max_time]
        self.data = self.data[min_time:max_time, :, min_x:max_x, min_y:max_y]
        self.lat = self.lat[min_x:max_x, min_y:max_y]
        self.lon = self.lon[min_x:max_x, min_y:max_y]


    def write_curv_grid(self, time, lon, lat, u, v, atts, ofn):
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
        nc_u   = nc.createVariable('water_u','f4',('time','y','x'))
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
        print "closing netcdf file"
        nc.close()
        
    def write_curv_grid_wind(self, time, lon, lat, u, v, atts, ofn):
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
        nc_u   = nc.createVariable('air_u','f4',('time','y','x'))
        nc_v = nc.createVariable('air_v','f4',('time','y','x'))
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
        print "closing netcdf file"
        nc.close()


## utilities:

def TD2float_days(td):
    return td.days + td.seconds / (24.0 * 3600)  

        
if __name__ == "__main__":
    
#    filename = "SampleData/SURF_GOM27_10"
#    grid = 'coarse'
#    filename = "SampleData/SURF_GOM27_82_02"
#    grid = 'fine'
#    filename = "SampleData/SURF_GOM27_50_02"
#    grid = 'fine'

    filename = 'C:/workingdisk/MMS_Data/Oey_atl_2003/SURF_ATL27_06.DAT'
    grid = 'coarse'

    reader = data_reader(filename, grid)

#    reader.read_all()
#    reader.to_netcdf_wind(filename+"_wind.nc")

    reader.read_times()

    reader.read_data()
    reader.read_grid()

#    raise Exception("I'm stopping now")

    
    print "start time:", reader.times[0]
    print "end time:", reader.times[-1]
    data = reader.data
    print "data shape:", data.shape
    print "some data is:", data[0,:,150,150]
    print "some stats:"
    print "min latitude:", reader.lat.min()
    print "max latitude:", reader.lat.max()
    print "min longitude:", reader.lon.min()
    print "max longitude:", reader.lon.max()
    for field, index in reader.fields.items():
        print "max %s"%field, data[:,index,:,:].max()
        print "min %s"%field, data[:,index,:,:].min()
    
    
