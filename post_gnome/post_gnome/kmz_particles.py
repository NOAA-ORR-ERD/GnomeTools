#!/usr/bin/env python

"""
Code for workign with particle fiels in mkz

Only handles reading for now
"""
# for py2/3 compatibility
from __future__ import absolute_import, division, print_function, unicode_literals

import nc_particles

from nc_particles import file_attributes, var_attributes

## variables used to support the stucture of the file, rather than data
## used to remove them from the list of available data variables
SPECIAL_VARIABLES = ['time','particle_count']


class Writer(object):
    def __init__ (self,
                  filename,
                  num_timesteps=None,
                  ref_time=None,
                  file_attributes=file_attributes,
                  var_attributes=var_attributes
                  ):

        """
        create a kmz_particle file Writer

        Creates the kml / kmz file, opens it for writing,
        writes the global attributes, creates required variables etc.

        :param filename: name of kmz file to open - if it exists,
                         it will be written over!

        :param num_timesteps=None: number of timesteps that will be output. Must be defined for netcdf3.
                                   Can be None for netcdf4
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
        else:
            # open a new one
            self.nc = netCDF4.Dataset(nc_file)

        time = self.nc.variables['time']
        units = time.getncattr('units')
        self.times = netCDF4.num2date(time[:], units)
        self.time_units = units

        self.particle_count = self.nc.variables['particle_count']
        # build the index:
        self.data_index = np.zeros((len(self.times) + 1,), dtype=np.int32)
        self.data_index[1:] = np.cumsum(self.particle_count)

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
        try:
            self.nc.close()
            print ("netcdf file closed")
        except RuntimeError:
            # just in case it isn't still open
            pass

    def __del__(self):
        """ make sure to close the netcdf file """
        self.close()