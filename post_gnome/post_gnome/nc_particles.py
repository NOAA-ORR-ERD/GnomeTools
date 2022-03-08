#!/usr/bin/env python

# for py2/3 compatibility
from __future__ import print_function, absolute_import, division,  unicode_literals


"""
Module for manipulating netcdf (and other) particle files

This is useful for working with particle file output from GNOME,
and is a test case for working with what hopefully will be a
CF standard (or SOME standard..)

"""

# from datetime import datetime

import numpy as np

import netCDF4

# default attributes -- can be updated by user later.

file_attributes = {'conventions': "CF-1.6",
                   'title': "Sample data/file for particle trajectory format",
                   'institution': "NOAA Emergency Response Division",
                   'source': "example data from nc_particles",
                   'history': "Evolved with discussion on CF-metadata listserve",
                   'references': '',
                   'comment': 'Some simple test data',
                   'CF:featureType': "particle_trajectory" ,
                   }

# variable attributes for the standard variables
# these will be used when the data are created.
# you can use a different set if you want.

var_attributes = {'time': {'long_name': 'time since the beginning of the simulation',
                           'standard_name': 'time',
                           'calendar': 'gregorian',
                           'standard_name': 'time',
                           'comment': 'unspecified time zone',
                           # units will get set based on data
                           },
                  'particle_count': {'units': '1',
                                     'long_name': 'number of particles in a given timestep',
                                     'ragged_row_count': 'particle count at nth timestep',
                                     },
                  'longitude': {'long_name': 'longitude of the particle',
                                'standard_name': 'longitude',
                                'units': 'degrees_east',
                                },
                  'latitude': {'long_name': 'latitude of the particle',
                               'standard_name': 'latitude',
                               'units': 'degrees_north',
                               },
                  'depth': {'long_name': 'particle depth below sea surface',
                            'standard_name': 'depth',
                            'units': 'meters',
                            'axis': 'z positive down',
                            },
                  'mass': {'long_name': 'mass of particle',
                           'units': 'grams',
                           },
                  'age': {'long_name': 'age of particle from time of release',
                          'units': 'seconds',
                          },
                  'status_code': {'long_name': 'particle status code',
                                  'valid_range': (0, 5),
                                  'flag_values': "7 12 0 10 2 3",
                                  'flag_meanings': "0: not_released, 2: in_water, 3: on_land, 7: off_maps, 10: evaporated, 12: to_be_removed,"
                                  },
                  'id': {'long_name': 'particle ID',
                         },
                  }

# alternate names:
var_attributes['lon'] = var_attributes['longitude']
var_attributes['lat'] = var_attributes['latitude']

# variables used to support the stucture of the file, rather than data
# used to remove them from the list of available data variables
SPECIAL_VARIABLES = ['time', 'particle_count']


class Writer(object):
    def __init__(self,
                 filename,
                 num_timesteps=None,
                 ref_time=None,
                 file_attributes=file_attributes,
                 var_attributes=var_attributes,
                 nc_version=4,
                 ):

        """
        create a nc_particle file Writer

        Creates the netccdf file, opens it for writing,
        writes the global attributes, creates required variables etc.

        :param filename: name of netcdf file to open - if it exists,
                         it will be written over!

        :param num_timesteps=None: number of timesteps that will be output. Must be
                                   defined for netcdf3. Can be None for netcdf4
        :type num_timesteps: integer

        :param ref_time=None: reference time for time units (i.e. seconds since..).
                              If None, the first time used will be used.
        :type ref_time: datetime object

        :param file_attributes: keys and values for teh file-level attributes.
                                Defaults to the set defined in this module.
        :type file_attributes: dict

        :param var_attributes: dist of variable names, and the keys and values for variable
                               attributes.
                               Defaults to the set defined in this module.
        :type var_attributes: dict

        :param nc_version=3: version of netcdf to use -- must be 3 or 4. If 4, some extra
                             features are enabled.
        :type nc_version: integer
        """

        self.num_timesteps = num_timesteps
        self.ref_time = ref_time

        self.file_attributes = file_attributes
        self.var_attributes = var_attributes
        self.nc = None

        try:
            nc_version = int(nc_version)
        except ValueError:
            raise ValueError("nc_format must be 3 or 4")
        if nc_version == 3:
            format = 'NETCDF3_CLASSIC'
        elif nc_version == 4:
            format = 'NETCDF4'
        else:
            raise ValueError("nc_format must be 3 or 4")

        if nc_version == 3 and self.num_timesteps is None:
            raise ValueError("You must specify num_timesteps when using netcdf3")

        nc = netCDF4.Dataset(filename, 'w', format=format)
        self.nc = nc

        # Global attributes
        for (name, value) in self.file_attributes.items():
            setattr(nc, name, value)

        # Dimensions
        nc.createDimension('time', self.num_timesteps)
        nc.createDimension('data', None)

        # required variables
        time = nc.createVariable('time', np.int32, ('time',))
        for name, value in self.var_attributes['time'].items():
            time.setncattr(name, value)
        # make sure there are some units there
        # this will get overwritten when the proper reference time is known
        if self.ref_time is None:
            # make sure there are some units there
            # this will get overwritten when the proper reference time is known
            time.units = "seconds since 2016-01-01T00:00:00"
        else:
            time.units = 'seconds since {0}'.format(self.ref_time.isoformat())

        pc = nc.createVariable('particle_count', np.int32, ('time',))
        for name, value in self.var_attributes['particle_count'].items():
            pc.setncattr(name, value)
        self.time_var = time

        self.num_data = 0
        self.current_timestep = 0

    def write_timestep(self, timestamp, data):
        """
        write the data for a timestep

        :param timestamp: the time stamp of the timestep
        :type timestamp: datetime object

        :param data: dict of data arrays -- all parameters for a single time step
        :type data: dict

        Note: it is assumed that the timesteps will be written sequentially,
              and that the variables will not change after the first timestep
              is written.
        """

        nc = self.nc
        if self.current_timestep == 0:
            # create the variables and add attributes
            # set the time units:
            if self.ref_time is None:
                self.ref_time = timestamp
            nc.variables['time'].units = 'seconds since {0}'.format(self.ref_time.isoformat())
            for key, val in data.iteritems():
                val = np.asarray(val)
                var = nc.createVariable(key, datatype=val.dtype, dimensions=('data'))
                # if it's a standard variable, add the attributes
                if key in self.var_attributes:
                    for name, value in self.var_attributes[key].items():
                        var.setncattr(name, value)

        particle_count = len(data.itervalues().next())  # length of an arbitrary array
        nc.variables['particle_count'][self.current_timestep] = particle_count
        nc.variables['time'][self.current_timestep] = (timestamp - self.ref_time).total_seconds()
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
        if self.nc is not None:
            try:
                self.nc.close()
            except RuntimeError:
                # just in case it isn't still open
                pass

    def __del__(self):
        """ make sure to close the netcdf file """
        self.close()


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
            self.manage_dataset = False
        else:
            # open a new one
            self.nc = netCDF4.Dataset(nc_file)
            self.manage_dataset = True

        time = self.nc.variables['time']
        units = time.getncattr('units')
        self.times = netCDF4.num2date(time[:], units)
        self.time_units = units

        self.particle_count = self.nc.variables['particle_count']
        # build the index:
        self.data_index = np.zeros((len(self.times) + 1,), dtype=np.int32)
        self.data_index[1:] = np.cumsum(self.particle_count)
        self.global_atttributes = {name: self.nc.getncattr(name) for name in self.nc.ncattrs()}

    @property
    def variables(self):
        """
        return the names of all the variables associated with the particles
        """
        return [var for var in self.nc.variables.keys() if var not in SPECIAL_VARIABLES]

    def __str__(self):
        return ("nc_particles Reader object:\n"
                "variables: {}\n"
                "number of timesteps: {}\n"
                ).format(self.variables, len(self.times))

    def get_all_timesteps(self, variables=['latitude', 'longitude']):
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
                ind2 = self.data_index[i + 1]
                data[var].append(self.nc.variables[var][ind1:ind2])
        return data

    def get_units(self, variable):
        """
        return the units of the given variable

        :param variable: name of the variable for which the units are required
        :type variable: string
        """
        return self.nc.variables[variable].units

    def get_attributes(self, variable):
        """
        return all the attributes of the given variable

        :param variable: name of the variable for which the attributes are required
        :type variable: string
        """
        var = self.nc.variables[variable]
        return {name: var.getncattr(name) for name in var.ncattrs()}

    def get_timestep(self, timestep, variables=['latitude', 'longitude']):
        """
        returns the requested variables data from a given timestep as a
        dictionary keyed by the variable names

        :param timestep: Index of the timestep you want the data for.
        :type timestep: int

        :param variables: The variables desired as a list string names.
                          Defaults to ['latitude','longitude']
        :type variables: list of strings

        :returns data: returns a dict of arrays -- the keys are the
                       variable names, and the values are numpy arrays
                       of the data.
        """
        ind1, ind2 = self.data_index[timestep:timestep + 2]
        return {var: self.nc.variables[var][ind1:ind2] for var in variables}

    def get_individual_trajectory(self, particle_id, variables=['latitude', 'longitude']):
        """
        returns the requested variables from trajectory of an individual particle

        :param particle_id: the id of the particle you want to track

        :param variables=['latitude', 'longitude']: which variables you want

        note: this is inefficient -- it has to read the entire file to get it.
        """
        indexes = np.where(self.nc.variables['id'][:] == particle_id)
        data = {}
        for var in variables:
            data[var] = self.nc.variables[var][indexes]
        return data

    def close(self):
        """
        close the netcdf file
        """
        if self.nc.isopen():
            self.nc.close()


    def __del__(self):
        """
        Make sure to close the netcdf file

        But only if this instance opened it in the first place
        """

        if self.manage_dataset:
            self.close()
