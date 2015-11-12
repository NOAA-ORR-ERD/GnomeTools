#!/usr/bin/env python

"""
module for manipulating netcdf particle files

This is useful for working with particle file output from GNOME,
and is a test case for working with what hopefully will be a
CF standard (or SOME standard..)

"""
from __future__ import division #Change the / operator to ensure true division throughout (Zelenke).
from datetime import datetime

import numpy as np

import netCDF4

file_attributes = {'conventions' : "CF-1.6",
                   'title' : "Sample data/file for particle trajectory format",
                   'institution' : "NOAA Emergency Response Division",
                   'source' : "example data from nc_particles",
                   'history' : "Evolved with discussion on CF-metadata listserve",
                   'references' : '',
                   'comment' : 'Some simple test data',
                   'CF:featureType' : "particle_trajectory" ,
                   }

## variable attributes for the standard variables
## these will be used when the data are created.
## you can use a different set if you want.

var_attributes = {'time': {'long_name':'time since the beginning of the simulation',
                           'standard_name':'time',
                           'calendar':'gregorian',
                           'standard_name':'time',
                           'comment':'unspecified time zone',
                           #units will get set based on data
                            },
                  'particle_count': {'units':'1',
                                     'long_name':'number of particles in a given timestep',
                                     'ragged_row_count':'particle count at nth timestep',
                                     },
                  'longitude': {'long_name':'longitude of the particle',
                                'standard_name':'longitude',
                                'units':'degrees_east',
                            },
                  'latitude': {'long_name':'latitude of the particle',
                           'standard_name':'latitude',
                            'units':'degrees_north',
                            },
                  'depth': {'long_name':'particle depth below sea surface',
                        'standard_name':'depth',
                        'units':'meters',
                        'axis':'z positive down',
                        },
                  'mass': {'long_name':'mass of particle',
                        'units':'grams',
                        },
                  'age': {'long_name':'age of particle from time of release',
                        'units':'seconds',
                        },
                  'status_code': {'long_name':'particle status code',
                                  'valid_range': (0, 5),
                                  'flag_values': (1, 2, 3, 4),
                                  'flag_meanings':'on_land off_maps evaporated below_surface'
                                  },
                  'id': {'long_name':'particle ID',
                        },
}
# alternate names:
var_attributes['lon'] = var_attributes['longitude']
var_attributes['lat'] = var_attributes['latitude']

## variables used to support the stucture of the file, rather than data
## used to remove them from the list of available data variables
SPECIAL_VARIABLES = ['time','particle_count']

class Writer(object):
    def __init__ (self, filename,
                        num_timesteps,
                        ref_time,
                        file_attributes=file_attributes,
                        var_attributes=var_attributes):

        """
        create a nc_particle file Writer

        Creates the netccdf file, opens it for writing,
        writes the global attributes, creates required variables etc.

        :param filename: name of netcdf file to open - if it exists,
                         it will be written over!

        :param num_timesteps: number of timesteps that will be output
        :type num_timesteps: integer

        :param ref_time: reference time for time units (i.e. seconds since..)
        :type ref_time: datetime object

        :param file_attributes: keys and values for teh file-level attributes.
                                Defaults to the set defined in this module.
        :type file_attributes: dict

        :param var_attributes: dist of variable names, and the keys and values for variable attributes.
                               Defaults to the set defined in this module.
        :type var_attributes: dict
        """

        self.num_timesteps = num_timesteps
        self.ref_time = ref_time

        self.file_attributes = file_attributes
        self.var_attributes = var_attributes

        nc = netCDF4.Dataset(filename, 'w', format='NETCDF3_CLASSIC')
        self.nc = nc

        ## Global attributes
        for (name, value) in self.file_attributes.items():
            setattr(nc, name, value)
        
        ## Dimensions
        nc.createDimension('time', self.num_timesteps )
        nc.createDimension('data', None)
        
        ## required variables
        time = nc.createVariable('time', np.int32, ('time',))
        time.units = 'seconds since {0}'.format(self.ref_time.isoformat())
        for name, value in var_attributes['time'].items():
            time.setncattr(name, value)

        pc = nc.createVariable('particle_count',np.int32, ('time',))
        for name, value in var_attributes['particle_count'].items():
            pc.setncattr(name, value)

        self.num_data = 0
        self.current_timestep = 0

    def write_timestep(self, timestep, data):
        """
        write the data for a timestep

        :param timestep: the time stamp of the timestep 
        :type timestep: datetime object

        :param data: dict of data arrays -- all parameters for a single time step
        :type data: dict

        Note: it is assumed that the timestpes will be written sequentially,
              and that the variables will not change after the first timestep
              is written.
        """

        nc = self.nc

        if self.current_timestep == 0:
            ## create the variables and add attributes
            for key, val in data.iteritems():
                var = nc.createVariable(key, datatype=val.dtype, dimensions=('data'))
                # if it's a standard variable, add the attributes
                if key in var_attributes:
                    for name, value in var_attributes[key].items():
                        var.setncattr(name, value)

        particle_count = len(data.itervalues().next()) # length of an arbitrary array
        nc.variables['particle_count'][self.current_timestep] = particle_count
        nc.variables['time'][self.current_timestep] = (timestep-self.ref_time).total_seconds()
        self.current_timestep += 1
        for key, val in data.iteritems():
            var = nc.variables[key]
            if len(val) != particle_count:
                raise ValueError("All data arrays must be the same length")
            var[self.num_data:] = val
        self.num_data += particle_count


    def close(self):
        """
        close the netcdf file
        """
        self.nc.close()

class Reader(object):
    """
    Class to handle reading a nc_particle file

    (such as those written by GNOME or the Writer class above)
    """
    def __init__(self, nc_file):
        """
        initialize a file reader.

        :param nc_file: the netcdf file to read. If a netCDF4 Dataset, it will be used,
                        if a string, a new netCDF Dataset will be opened for reading
                        using that filename
        :type nc_file: string or netCDF4 Dataset object

        """

        if type(nc_file) == netCDF4.Dataset:
            # already open -- just use it
            self.nc = nc_file
        else:
            # open a new one
            self.nc = netCDF4.Dataset(nc_file)

        time = self.nc.variables['time']
        units = time.getncattr('units')
        self.times = netCDF4.num2date(time[:], units)
        self.time_units = units
                
        self.particle_count = self.nc.variables['particle_count']
        # build the index:
        self.data_index = np.zeros((len(self.times)+1,), dtype=np.int32 )
        self.data_index[1:] = np.cumsum(self.particle_count)

    @property
    def variables(self):
        """
        return the names of all the variables associated with the particles
        """
        return [ var for var in self.nc.variables.keys() if var not in SPECIAL_VARIABLES]


    def get_all_timesteps(self, variables=['latitude','longitude']):
         """
         returns the requested variables data from all timesteps as a
         dictionary keyed by the variable names

         :param variables: the variables desired as a list string names.
                           Defaults to ['latitude','longitude']
         :type variables: list of strings

         :returns data: returns a dict of arrays -- the keys are the
                        variable names, and the values are numpy arrays
                        of the data. The arrays are the flattened ragged
                        array of data.
         """
         data = {}
         for var in variables:
            data[var] = []
            for i in range(len(self.times)):
                ind1 = self.data_index[i]
                ind2 = self.data_index[i+1]
                data[var].append( self.nc.variables[var][ind1:ind2] )      
         return data
   
    def get_timestep_single_var(self, timestep, variable):
        ind1 = self.data_index[timestep]
        ind2 = self.data_index[timestep+1]
        return self.nc.variables[variable][ind1:ind2]
   
    def get_units(self, variable):
        """
        return the units of the given variable

        :param variable: name of the variable for which the units are required
        :type variable: string
        """
        return self.nc.variables[variable].units

    def get_timestep(self, timestep, variables=['latitude','longitude']):
        """
        returns the requested variables data from a given timestep as a
        dictionary keyed by the variable names

        :param variables: The variables desired as a list string names.
                          Defaults to ['latitude','longitude']
        :type variables: list of strings

        :returns data: returns a dict of arrays -- the keys are the
                       variable names, and the values are numpy arrays
                       of the data.
        """
        ind1, ind2 = self.data_index[timestep:timestep+2]
        return {var:self.nc.variables[var][ind1:ind2] for var in variables}
        
    def get_individual_trajectory(self, particle_id, variables=['latitude','longitude']):
        """
        returns the requested variables from trajectory of an individual particle
        
        note: this is inefficient -- it has to read the entire file to get it.
        """
        indexes = np.where(self.nc.variables['id'][:] == particle_id)
        data = {}
        for var in variables:
            data[var] = self.nc.variables[var][indexes]      
        return data
    
    def close(self):
        self.nc.close()
        
        
    
                                             

