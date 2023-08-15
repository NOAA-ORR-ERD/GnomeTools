# -*- coding: utf-8 -*-
"""
Assorted utilities for plotting grid_plume results

@author: brian.zelenke
"""

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import shutil
import os
import datetime
from dateutil import parser

#==============================================================================
# Section #1:  "Convenience" functions.
#==============================================================================


def load_gridplume_netcdf(path2file):
    '''
    Load the data from a NetCDF file output by grid_plume.py into separate
    variables.
    Input:
        "path2file" is the full path to a NetCDF file output by grid_plume.py.
            The NetCDF file must contain the following variables in the units given:
                "lat" in units of decimal degrees north;
                "lon" in units of decimal degrees east;
                "depth" in units of meters;
                "time" in units of seconds since a date-time (e.g., Unix time)
                "concentration" in units of kilograms per cubic meter.
    Output:
        "lat" in units of decimal degrees north;
        "lon" in units of decimal degrees east;
        "depth" in units of meters;
        "time" in units of seconds since a date-time
        "data" is the open NetCDF dataset.  The concentration variable can be
            read from this by the plotting function, either piecemeal or in
            its entirety as system memory allows.
        "dt" is a Python datetime object of the times in "time"
    Example:
        lat, lon, depth, time, data, dt = plot_results.load_gridplume_netcdf( 'C:/Users/brian.zelenke/Documents/GNOME/UXO/GlenWatabayashi/UXO_ROMS/Zelenke_test/20120315_GridPlumeOutput3Days10kg1000les_uh.nc' )
    '''

    data = nc.Dataset(path2file, 'r')

    if data.variables['concentration'].getncattr('units').lower() != 'kg/m^3':
        raise ValueError(
            'Concentration must be in units of kg/m^3. Instead, units were %s.'
            % data.variables['concentration'].getncattr('units'))

    if data.variables['lat'].getncattr('units').lower() != 'degrees_north':
        raise ValueError(
            'Latitude must be in units of degrees_north. Instead, units were %s.'
            % data.variables['lat'].getncattr('units'))

    if data.variables['lon'].getncattr('units').lower() != 'degrees_east':
        raise ValueError(
            'Longitude must be in units of degrees_east. Instead, units were %s.'
            % data.variables['lon'].getncattr('units'))

    if data.variables['depth'].getncattr('units').lower() != 'meters':
        raise ValueError(
            'Depth must be in units of meters. Instead, units were %s.' %
            data.variables['depth'].getncattr('units'))

    if data.variables['time'].getncattr(
            'units').lower()[0:13] != 'seconds since':  #Is this necessary?
        raise ValueError(
            'Time must be in units of seconds since a date-time (e.g., Unix time). Instead, units were %s.'
            % data.variables['time'].getncattr('units'))

    lon = data.variables['lon'][:]
    lat = data.variables['lat'][:]
    depth = data.variables['depth'][:]
    time = data.variables['time'][:]

    #Create datetime object from elapsed seconds stored in NetCDF.
    EpochString = data.variables['time'].getncattr('units').lower()[-19:]
    EpochDateTime = parser.parse(EpochString)
    dt = [
        EpochDateTime + datetime.timedelta(seconds=int(time[i]))
        for i in range(len(time))
    ]

    #Concentration isn't returned.
    #The following line is commented-out due to the concentration matrix
    #frequently surpassing system memory.

    #concentration = data.variables['concentration'][:]

    #Nor is a NetCDF variable of concentration returned.
    #The following line is commented-out because if the NetCDF dataset this
    #points to is closed (which it should be after use) the NetCDF ID this
    #variable refers to is lost.  Instead the whole NetCDF dataset is returned.
    #Keeping the dataset around is also handy for copying attributes into
    #derivative NetCDF files of output.

    #NetCDFconcentration = data.variables['concentration']

    return (
        lat,
        lon,
        depth,
        time,
        data,
        dt,
    )


def convert_concentration(array, units):
    '''
    Convert the input array of concentrations from SI units of kilograms per
    cubic meter to desired units.
    Input:
        "array" is an array of concentration values in units of kilograms per
            cubic meter.
        "units" are the units concentration will be returned in; either:
            "ppt" for parts per thousand (a.k.a. kilograms per cubic meter)
            "%" for percent
            "fraction" for fraction (decimal)
            "micro-g/l" for microgram per liter
            "mg/l" for milligrams per liter
            "ng/l" for nanograms per liter
            "ppb" for parts per billion
            "ppm" for parts per million
            "pptrillion" for parts per trillion
            "lb/ft^3" for pounds per cubic foot
    Output:
        "concentrations" is an array of concentration values in units per "units".
    Example:
        concentrations=plot_results.convert_concentration(concentration,'ppb')
    '''

    #Specify factor for unit conversion.
    #FixMe:  From hazpy import unit_conversion (see what's within unit_conversion... so could do unit_conversion.[something])
    if units.lower() == 'ppt':
        factor = 1.0
    elif units.lower() == '%':
        factor = 0.1
    elif units.lower() == 'fraction':
        factor = 0.001
    elif units.lower() == 'micro-g/l' or units.lower() == 'ppb':
        factor = 1.e6
    elif units.lower() == 'mg/l' or units.lower() == 'ppm':
        factor = 1.e3
    elif units.lower() == 'ng/l' or units.lower() == 'pptrillion':
        factor = 1.e9
    elif units.lower() == 'lb/ft^3':
        factor = 0.0624
    else:
        raise ValueError(
            'Concentration units "%s" are not supported.\nConcentration can be output in units of "ppt", "%%", "fraction", "micro-g/l", "mg/l", "ng/l", "ppb", "ppm", "pptrillion", or "lb/ft^3".'
            % units)

    #Convert concentrations from one unit to another.
    concentrations = np.multiply(array, factor)
    return concentrations


def find_nearest(array, value):
    '''
    Find the number in the array nearest to the value given.
    Input:
        "array" is a NumPy array of numbers.
        "value" is a number.
    Output:
        "array" is the number in the array closest to the given value.  If more
            than one number in the array is equally close to the value given,
            the first of these numbers to appear in the array will be returned.
    Example:
        num=plot_results.find_nearest(np.array([7,2,5]),6)
        num=7
    '''
    idx = (np.abs(array - value)).argmin()
    return array[idx]


#==============================================================================
# Section #2:  Calculation functions.
#==============================================================================


def max_conc_grids(path2file, units, path4output=None, levels=None):
    '''
    Create a grid of maximum concentrations measured at each grid-point during the sampling time, irrespective of depth.
    Input:
        "path2file" is the full path to a NetCDF file output by grid_plume.py.
            The NetCDF file must contain the following variables in the units given:
                "lat" in units of decimal degrees north;
                "lon" in units of decimal degrees east;
                "depth" in units of meters;
                "time" in units of seconds since a date-time (e.g., Unix time)
                "concentration" in units of kilograms per cubic meter.
         "units" are the units concentration will be plotted in; either:
            "ppt" for parts per thousand (a.k.a. kilograms per cubic meter)
            "%" for percent
            "fraction" for fraction (decimal)
            "micro-g/l" for microgram per liter
            "mg/l" for milligrams per liter
            "ng/l" for nanograms per liter
            "ppb" for parts per billion
            "ppm" for parts per million
            "pptrillion" for parts per trillion
            "lb/ft^3" for pounds per cubic foot
        "path4output" is optional and is the full path to the file where you
            want to save the results as either a comma-separated-value (CSV)
            spreadsheet or NetCDF-formatted file. If "path4output" ends in
            ".csv" a CSV-formatted spreadsheet will be output.  If
            "path4output" ends in ".nc" a NetCDF file will be returned.
    Output:
        "Lon2D" is a two-dimensional rectangular grid of the longitudes from
            "path2file".
        "Lat2D" is a two-dimensional rectangular grid of the latitudes from
            "path2file".
        "MaxConc2D" is a two-dimensional rectangular grid of the maximum:
            concentration found within the water column during the sampling
            time at each Lat2D/Lon2D grid-point.  If no concentration value was
            recorded at a Lat2D/Lon2D grid-point in the time-series,
            "MaxConc2D" has a value of NaN at that coordinate.
        Optionally either a CSV or NetCDF file of the results, depending on the
            file extension specified in "path4output".
    Example:
        Lon2D, Lat2D, MaxConc2D = plot_results.max_conc_grids( 'C:/Users/brian.zelenke/Documents/GNOME/Tests/GNOME_3d_test/2013-03-14_PresentationToDebbiePayton/20130312_GridPlumeOutput_Zelenke_All_spills_500m.nc', 'ppm' )

    '''

    lat, lon, _, time, data, _ = load_gridplume_netcdf(path2file)

    # Create 2-D rectangular grids for plotting.
    Lon2D, Lat2D = np.meshgrid(lon, lat)
    # Create 3-D rectangular grids for plotting (would first need to 'import scitools').
    # Lon3D,Lat3D,Depth3D=scitools.numpyutils.meshgrid(lon,lat,depth)

    # Evaluate concentration grid; even if it's too big to all fit into memory.
    try:
        if levels is None:

            MaxConcGrid = data.variables['concentration'][:]
        else:
            l_slice = slice(*levels)
            print("using:", l_slice)
            MaxConcGrid = data.variables['concentration'][:, l_slice, :, :]
        #Calculate the maximum concentration measured at each grid-point during the sampling time, irrespective of depth.
        MaxConcGrid = np.nanmax(
            MaxConcGrid, axis=0
        )  #Calculate the maximum during the time-series at each grid-point, for each time.
        MaxConcGrid = np.nanmax(
            MaxConcGrid, axis=0
        )  #Calculate the maximum during the time-series at each grid-point, across the depth bins.
    except ValueError:
        raise
        MaxConcGrid = np.zeros_like(Lat2D)
        for t_index in range(len(time)):
            print("Processing time-step %g of %g (%g%%)..." %
                  (t_index, len(time) - 1, 100 * (t_index / (len(time) - 1))))
            tMaxConc = data.variables['concentration'][
                t_index, :]  #Get concentration grid for a single time-step (hopefully the results aren't so big this [z,y,x] matrix won't fit into memory... if they are, could do another loop).
            tMaxConc = np.nanmax(
                tMaxConc, axis=0
            )  #Calculate the maximum during that time-step at each grid-point, across the depth bins.
            MaxConcGrid = np.fmax(MaxConcGrid, tMaxConc)

    #Mask-out areas with no concentration values.
    MaxConc2D = convert_concentration(
        MaxConcGrid,
        units)  #Convert concentration units (i.e., from units of kg/[m^3]).
    MaxConc2D[MaxConc2D == 0] = np.NaN

    if path4output != None:

        if (path4output.lower()[-3:]
                == ".nc") or (path4output.lower()[-4:]
                              == ".cdf") or (path4output.lower()[-7:]
                                             == ".netcdf"):

            #Initialize NetCDF file.
            nc_MaxConc2D = nc.Dataset(path4output, 'w', format="NETCDF4")

            #Specify NetCDF dimensions.
            nc_MaxConc2D.createDimension('nLatitudes', len(lat))  #Y
            nc_MaxConc2D.createDimension('nLongitudes', len(lon))  #X

            #Mesh coordinates.
            lats = nc_MaxConc2D.createVariable('latitude',
                                               np.double,
                                               ('nLatitudes', 'nLongitudes'),
                                               zlib=True)
            lats.standard_name = data.variables['lat'].getncattr(
                'standard_name')
            lats.long_name = "Latitude coordinate matrix."
            lats.units = data.variables['lat'].getncattr('units')
            lats[:] = Lat2D

            lons = nc_MaxConc2D.createVariable('longitude',
                                               np.double,
                                               ('nLatitudes', 'nLongitudes'),
                                               zlib=True)
            lons.standard_name = data.variables['lon'].getncattr(
                'standard_name')
            lons.long_name = "Longitude coordinate matrix."
            lons.units = data.variables['lon'].getncattr('units')
            lons[:] = Lon2D

            #Maximum concentrations at grid-points.
            concs = nc_MaxConc2D.createVariable('concentration',
                                                np.double,
                                                ('nLatitudes', 'nLongitudes'),
                                                zlib=True)
            concs.short_name = "mass_concentration_of_pollutant_in_sea_water"  #Note:  not a standard name.
            concs.long_name = "Matrix of maximum concentrations modeled at each grid coordinate during the sampling time, irrespective of depth."
            concs.units = units
            concs[:] = MaxConc2D

            #Set file attributes.
            #setattr( nc_MaxConc2D, "Conventions", "CF-1.6" )
            setattr(nc_MaxConc2D, "Title",
                    "GNOME GridPlume: Maximum Concentration")
            setattr(nc_MaxConc2D, "Institution", "USA/DOC/NOAA/NOS/ERD/TSSB")
            setattr(nc_MaxConc2D, "References", "brian.zelenke@noaa.gov")

            #Close writing to NetCDF file.
            nc_MaxConc2D.close()

        else:
            HeaderRow = "Latitude (\N{DEGREE SIGN}), Longitude (\N{DEGREE SIGN}), Maximum Concentration (%s)" % units
            database = np.column_stack(
                (Lat2D.flatten(), Lon2D.flatten(), MaxConc2D.flatten()))
            #database = np.ma.array( database, mask=np.isnan(database) ) #TODO:  Figure out how to make NaN values blank in CSV output (instead of "nan").
            np.savetxt(path4output,
                       database,
                       fmt='%g',
                       delimiter=",",
                       header=HeaderRow,
                       comments="")

    data.close()

    return (
        Lon2D,
        Lat2D,
        MaxConc2D,
    )


def point_conc(path2file, units, PointLat, PointLon, path4output=None):
    '''
    Calculate the time-series of the highest concentration modeled at any depth
    at the grid location closest to the coordinate given.
    Input:
        "path2file" is the full path to a NetCDF file output by grid_plume.py.
            The NetCDF file must contain the following variables in the units given:
                "lat" in units of decimal degrees north;
                "lon" in units of decimal degrees east;
                "depth" in units of meters;
                "time" in units of seconds since a date-time (e.g., Unix time)
                "concentration" in units of kilograms per cubic meter.
        "units" are the units concentration will be plotted in; either:
            "ppt" for parts per thousand (a.k.a. kilograms per cubic meter)
            "%" for percent
            "fraction" for fraction (decimal)
            "micro-g/l" for microgram per liter
            "mg/l" for milligrams per liter
            "ng/l" for nanograms per liter
            "ppb" for parts per billion
            "ppm" for parts per million
            "pptrillion" for parts per trillion
            "lb/ft^3" for pounds per cubic foot
        "PointLat" is the latitude coordinate in decimal degrees (+=northern
            hemisphere, -=southern hemisphere) of your location of interest
            that falls within the area modeled with GridPlume.
        "PointLon" is the longitude coordinate in decimal degrees (+=eastern
            hemisphere, -=western hemisphere) of your location of interest that
            falls within the area modeled with GridPlume.
        "path4output" is optional and is the full path to the file where you
            want to save the results as either a comma-separated-value (CSV)
            spreadsheet or NetCDF-formatted file. If "path4output" ends in
            ".csv" a CSV-formatted spreadsheet will be output.  If
            "path4output" ends in ".nc" a NetCDF file will be returned.
    Output:
        "MaxConcAtPt" is the time-series of the highest concentration modeled
            at any depth at the grid location closest to the coordinate given,
            returned in the "units" selected.
        "time" is the aforementioned "time" variable from "path2file".
        "ptLat" is the latitude, in units of decimal degrees, of the grid point
            closest to the input coordinate.
        "ptLon" is the longitude, in units of decimal degrees, of the grid
            point closest to the input coordinate.
        "DepthOfMaxConcAtPt" is the time-series of depths at which
            "MaxConcAtPt" occurred.  If the maximum concentration at a
            time-step occured at more than one depth, the shallowest will be
            returned within "DepthOfMaxConcAtPt".
        Optionally either a CSV or NetCDF file of the results, depending on the
            file extension specified in "path4output".
    Example:
        MaxConcAtPt, time, ptLat, ptLon, DepthOfMaxConcAtPt = plot_results.point_conc( 'C:/Users/brian.zelenke/Documents/GNOME/UXO/GlenWatabayashi/UXO_ROMS/Zelenke_test/20120315_GridPlumeOutput3Days10kg1000les_uh.nc', 'ppb', 21.453471, -158.204816 )
    '''

    #Load the NetCDF file output by grid_plume.py.
    lat, lon, depth, time, data, dt = load_gridplume_netcdf(path2file)

    #Create 2-D rectangular grids for plotting.
    Lon2D, Lat2D = np.meshgrid(lon, lat)

    #Find the indicies of the coordinate within the 2-D grid closest to the point given.
    ptLon = find_nearest(lon, PointLon)
    ptLat = find_nearest(lat, PointLat)
    ij = np.argwhere(np.logical_and(Lon2D == ptLon, Lat2D == ptLat) == True)

    #Calculate the maximum concentration at each time-step across the depth bins at the given point in the grid.
    ConcAtPt = data.variables['concentration'][:, :, ij[0, 0], ij[0, 1]]
    MaxConcAtPt = np.nanmax(ConcAtPt, axis=1)

    zi = np.nanargmax(
        ConcAtPt, axis=1
    )  #Depth indices of grid-point's maximum concentrations at each time-step.
    DepthOfMaxConcAtPt = depth[
        zi]  #Depths of maximum concentrations at grid-point.

    MaxConcAtPt = convert_concentration(
        MaxConcAtPt,
        units)  #Convert concentration units (i.e., from units of kg/[m^3]).

    if path4output != None:

        if (path4output.lower()[-3:]
                == ".nc") or (path4output.lower()[-4:]
                              == ".cdf") or (path4output.lower()[-7:]
                                             == ".netcdf"):

            #Initialize NetCDF file.
            nc_PointConc2D = nc.Dataset(path4output, 'w', format="NETCDF4")

            #Specify NetCDF dimensions.
            nc_PointConc2D.createDimension('nTimes', len(time))
            nc_PointConc2D.createDimension('One', 1)

            #T, X, Y, and Z coordinates.
            times = nc_PointConc2D.createVariable('time',
                                                  np.double, ('nTimes', ),
                                                  zlib=True)
            times.standard_name = data.variables['time'].getncattr(
                'standard_name')
            times.long_name = data.variables['time'].getncattr('long_name')
            times.units = data.variables['time'].getncattr('units')
            times[:] = time

            lats = nc_PointConc2D.createVariable('latitude',
                                                 np.double, ('One', ),
                                                 zlib=True)
            lats.standard_name = data.variables['lat'].getncattr(
                'standard_name')
            lats.long_name = data.variables['lat'].getncattr('long_name')
            lats.units = data.variables['lat'].getncattr('units')
            lats[:] = ptLat

            lons = nc_PointConc2D.createVariable('longitude',
                                                 np.double, ('One', ),
                                                 zlib=True)
            lons.standard_name = data.variables['lon'].getncattr(
                'standard_name')
            lons.long_name = data.variables['lon'].getncattr('long_name')
            lons.units = data.variables['lon'].getncattr('units')
            lons[:] = ptLon

            depths = nc_PointConc2D.createVariable('depth',
                                                   np.double, ('nTimes', ),
                                                   zlib=True)
            depths.standard_name = data.variables['depth'].getncattr(
                'standard_name')
            depths.long_name = "Depths where maximum concentration at coordinate was modeled to occur."
            depths.units = data.variables['depth'].getncattr('units')
            depths[:] = DepthOfMaxConcAtPt

            #Maximum concentrations at grid point.
            concs = nc_PointConc2D.createVariable('concentration',
                                                  np.double, ('nTimes', ),
                                                  zlib=True)
            concs.short_name = "mass_concentration_of_pollutant_in_sea_water"  #Note:  not a standard name.
            concs.long_name = "Time-series of maximum concentration modeled at coordinate."
            concs.units = units
            concs[:] = MaxConcAtPt

            #Set file attributes.
            #setattr( nc_MaxConc2D, "Conventions", "CF-1.6" )
            setattr(nc_PointConc2D, "Title",
                    "GNOME GridPlume: Maximum Concentrations at Coordinate")
            setattr(nc_PointConc2D, "Institution", "USA/DOC/NOAA/NOS/ERD/TSSB")
            setattr(nc_PointConc2D, "References", "brian.zelenke@noaa.gov")

            #Close writing to NetCDF file.
            nc_PointConc2D.close()

        else:
            HeaderRow = "Time, Maximum Concentration (%s), Depth of Maximum (m), Latitude for All Values (\N{DEGREE SIGN}), Longitude for All Values (\N{DEGREE SIGN})\n%s,%g,%g,%g,%g" % (
                units, dt[0], MaxConcAtPt[0], DepthOfMaxConcAtPt[0], ptLat,
                ptLon)
            database = np.column_stack(
                (dt[1:], MaxConcAtPt[1:], DepthOfMaxConcAtPt[1:]))
            np.savetxt(path4output,
                       database,
                       fmt=['%s', '%g', '%g'],
                       delimiter=",",
                       header=HeaderRow,
                       comments="")

    data.close()

    return (
        MaxConcAtPt,
        time,
        ptLat,
        ptLon,
        DepthOfMaxConcAtPt,
    )


def max_conc_in_plume(path2file, units, path4output=None, levels=None):
    '''
    Calculate the time-series of the highest concentration modeled at any depth
    anywhere within the given grid.
    Input:
        "path2file" is the full path to a NetCDF file output by grid_plume.py.
            The NetCDF file must contain the following variables in the units given:
                "lat" in units of decimal degrees north;
                "lon" in units of decimal degrees east;
                "depth" in units of meters;
                "time" in units of seconds since a date-time (e.g., Unix time)
                "concentration" in units of kilograms per cubic meter.
        "units" are the units concentration will be plotted in; either:
            "ppt" for parts per thousand (a.k.a. kilograms per cubic meter)
            "%" for percent
            "fraction" for fraction (decimal)
            "micro-g/l" for microgram per liter
            "mg/l" for milligrams per liter
            "ng/l" for nanograms per liter
            "ppb" for parts per billion
            "ppm" for parts per million
            "pptrillion" for parts per trillion
            "lb/ft^3" for pounds per cubic foot
        "path4output" is optional and is the full path to the file where you
            want to save the results as either a comma-separated-value (CSV)
            spreadsheet or NetCDF-formatted file. If "path4output" ends in
            ".csv" a CSV-formatted spreadsheet will be output.  If
            "path4output" ends in ".nc" a NetCDF file will be returned.
        :param levels=None: If you want a subset of the levels: tuple of indexes:
                            (0, 10) will get you the top 10
    Output:
        "MaxConcSeries" is the time-series of the highest concentration modeled
            at any depth anywhere within the given grid, returned in the
            "units" selected.  Should the concentration plume dissipate
            below detectable levels (viz. maximum concentration at zero), NaNs
            will be returned at remaining time-steps.
        "time" is the aforementioned "time" variable from "path2file".
        Optionally either a CSV or NetCDF file of the results, depending on the
            file extension specified in "path4output".
    Example:
        MaxConcSeries, time = plot_results.max_conc_in_plume( 'C:/Users/brian.zelenke/Documents/GNOME/UXO/GlenWatabayashi/UXO_ROMS/Zelenke_test/20120315_GridPlumeOutput3Days10kg1000les_uh.nc', 'ppb' )
    '''

    #Load the NetCDF file output by grid_plume.py.
    lat, lon, _, time, data, dt = load_gridplume_netcdf(path2file)

    #Calculate the maximum concentration measured during each time-step, irrespective of location or depth.
    try:
        MaxConcSeries = data.variables['concentration'][:]
        MaxConcSeries = np.nanmax(
            MaxConcSeries,
            axis=1)  #Calculate the maximum across the depth dimension.
        MaxConcSeries = np.nanmax(
            MaxConcSeries,
            axis=1)  #Calculate the maximum across the latitude dimension.
        MaxConcSeries = np.nanmax(
            MaxConcSeries,
            axis=1)  #Calculate the maximum across the longitude dimension.
    except ValueError:
        MaxConcSeries = np.zeros_like(time)  #Preallocate
        for t_index in range(len(time)):
            print("Processing time-step %g of %g (%g%%)..." %
                  (t_index, len(time) - 1, 100 * (t_index / (len(time) - 1))))
            tMaxConc = data.variables['concentration'][
                t_index, :]  #Get concentration grid for a single time-step (hopefully the results aren't so big this [z,y,x] matrix won't fit into memory... if they are, could do another loop).
            MaxConcSeries[t_index] = np.nanmax(tMaxConc)

    MaxConcSeries[
        MaxConcSeries ==
        0] = np.NaN  #Now just a time-series of maximums in the time dimension remains (less any times where the plume dissipated entirely.)
    MaxConcSeries = convert_concentration(
        MaxConcSeries,
        units)  #Convert concentration units (i.e., from units of kg/[m^3]).

    if path4output != None:

        if (path4output.lower()[-3:]
                == ".nc") or (path4output.lower()[-4:]
                              == ".cdf") or (path4output.lower()[-7:]
                                             == ".netcdf"):

            #Initialize NetCDF file.
            nc_PointConc2D = nc.Dataset(path4output, 'w', format="NETCDF4")

            #Specify NetCDF dimensions.
            nc_PointConc2D.createDimension('nTimes', len(time))

            #Time-steps.
            times = nc_PointConc2D.createVariable('time',
                                                  np.double, ('nTimes', ),
                                                  zlib=True)
            times.standard_name = data.variables['time'].getncattr(
                'standard_name')
            times.long_name = data.variables['time'].getncattr('long_name')
            times.units = data.variables['time'].getncattr('units')
            times[:] = time

            #Maximum concentrations modeled in grid.
            concs = nc_PointConc2D.createVariable('concentration',
                                                  np.double, ('nTimes', ),
                                                  zlib=True)
            concs.short_name = "mass_concentration_of_pollutant_in_sea_water"  #Note:  not a standard name.
            concs.long_name = "Time-series of maximum concentrations within model domain."
            concs.units = units
            concs[:] = MaxConcSeries

            #Set file attributes.
            #setattr( nc_MaxConc2D, "Conventions", "CF-1.6" )
            setattr(nc_PointConc2D, "Title",
                    "GNOME GridPlume: Maximum Concentrations Modeled")
            setattr(nc_PointConc2D, "Institution", "USA/DOC/NOAA/NOS/ERD/TSSB")
            setattr(nc_PointConc2D, "References", "brian.zelenke@noaa.gov")

            #Close writing to NetCDF file.
            nc_PointConc2D.close()

        else:
            HeaderRow = "Time, Maximum Concentration (%s)" % units
            database = np.column_stack((dt, MaxConcSeries))
            np.savetxt(path4output,
                       database,
                       fmt=['%s', '%g'],
                       delimiter=",",
                       header=HeaderRow,
                       comments="")

    data.close()

    return (
        MaxConcSeries,
        time,
    )


#==============================================================================
# Section #3:  Figure production functions.
#==============================================================================


def make_3ColorFig(Lon2D,
                   Lat2D,
                   masked_array,
                   units,
                   levels,
                   FontPtSize=14,
                   isGoogleEarthKMZ=False):
    '''
    Create a figure with concentration plotted as just three contours of "high-medium-low" colors.
    Input:
        "Lon2D" is a two-dimensional rectangular grid of the longitudes
            corresponding to the concentration values in "masked_array".
        "Lat2D" is a two-dimensional rectangular grid of the latitudes
            corresponding to the concentration values in "masked_array".
        "masked_array" is a masked array of a rectangular grid of
            concentrations with zero values masked-out.
        "units" are the units of the concentration values in masked_array; either:
            "ppt" for parts per thousand (a.k.a. kilograms per cubic meter)
            "%" for percent
            "fraction" for fraction (decimal)
            "micro-g/l" for microgram per liter
            "mg/l" for milligrams per liter
            "ng/l" for nanograms per liter
            "ppb" for parts per billion
            "ppm" for parts per million
            "pptrillion" for parts per trillion
            "lb/ft^3" for pounds per cubic foot
            This "units" string is just used only to label the colorbar.
        "levels" is a Numpy array of four values specifying (1) the bottom of
            the "low" end of the color scale, (2) the concentration where the
            color scale transitions from low to medium, (3) the concentration
            where the color scale transitions from medium to high, and (4) the
            top of the "high" end of the color scale.  Values must be in the
            same units as the concentrations in "masked_array".
        "FontPtSize" is optional and is the size of the font in points to use
            in labeling the plot.
        "isGoogleEarthKMZ" is optional and, if set to true, will return a
            figure without axes, ticks, or labels for use as an image overlay
            in Google Earth.
    Output:
        "csfig" is an instance of a contour figure that can be used to generate
            or save the graphic using MatPlotLib.
        "minX" is the lower limit of the x (longitude) axis.
        "maxX" is the upper limit of the x (longitude) axis.
        "minY" is the lower limit of the y (latitude) axis.
        "maxY" is the upper limit of the y (latitude) axis.
        "cbarfig" is a seperate figure instance that just contains the colorbar
            that goes with "csfig".  A figure instance is returned if
            "isGoogleEarthKMZ" is True; otherwise a NaN is returned.
    Example:
        csfig, minX, maxX, minY, maxY, cbarfig = plot_results.make_3ColorFig( Lon2D, Lat2D, masked_array, units, levels, FontPtSize=14, isGoogleEarthKMZ=False )
    '''

    levels = np.sort(levels[:])  #Quality control.
    if len(levels) != 4:
        raise ValueError('MAKE_3COLORFIG only supports three contour levels.')

    #Calculate the maximum and minimum concentrations.
    maxConc = np.nanmax(masked_array.data[:])
    minConc = np.nanmin(masked_array.data[:])

    #Get the extrema of the contour levels.
    MinConc2Plot = np.min(levels)
    MaxConc2Plot = np.max(levels)

    #Create the figure instance and axis.
    if isGoogleEarthKMZ:
        csfig = plt.figure(frameon=False)
        csax = csfig.add_axes((0, 0, 1, 1))
    else:
        csfig = plt.figure()
        csax = csfig.add_subplot(1, 1, 1)

    #Plot the contour on the axis.  Make extrema of colorbar arrow shaped if
    #concentration values exist beyond the limits specified.
    if MaxConc2Plot < maxConc and MinConc2Plot > minConc:
        cssurf = csax.contourf(Lon2D,
                               Lat2D,
                               masked_array,
                               levels,
                               colors=['blue', 'black', 'green'],
                               extend="both")
    elif MaxConc2Plot < maxConc and MinConc2Plot <= minConc:
        cssurf = csax.contourf(Lon2D,
                               Lat2D,
                               masked_array,
                               levels,
                               colors=['blue', 'black', 'green'],
                               extend="max")
    elif MaxConc2Plot >= maxConc and MinConc2Plot > minConc:
        cssurf = csax.contourf(Lon2D,
                               Lat2D,
                               masked_array,
                               levels,
                               colors=['blue', 'black', 'green'],
                               extend="min")
    else:
        cssurf = csax.contourf(Lon2D,
                               Lat2D,
                               masked_array,
                               levels,
                               colors=['blue', 'black', 'green'],
                               extend="neither")

    if isGoogleEarthKMZ:
        minX, maxX, minY, maxY = csax.axis()
        csax.set_axis_off()
        #csfig.savefig(path4image, transparent=True, bbox_inches=0, pad_inches=0)

        #Create a seperate figure instance that just has the colorbar that goes with the contourf plot.
        cbarfig = plt.figure()
        cscbar = cbarfig.colorbar(cssurf)
        cscbar.set_label('Concentration (%s)' % units,
                         fontsize=FontPtSize,
                         fontweight='bold')
        cscbar.ax.tick_params(
            labelsize=FontPtSize)  #Do you want to add a tick size parameter?
        cbarfig.delaxes(cbarfig.gca())
        #cbarfig.savefig(path4cbar, transparent=False, bbox_inches='tight')
    else:
        #cssurf.set_clim(MinConc2Plot,MaxConc2Plot) #Here "levels" sets the clim.
        #csax.set_aspect('equal') #Approximate a proper map projection via equal scaling of the plot's axes.
        csax.ticklabel_format(style='plain', useOffset=False, axis='both')
        csax.tick_params(labelsize=FontPtSize)
        csax.set_xlabel('Longitude ($\degree$)',
                        fontsize=FontPtSize,
                        fontweight='bold')
        csax.set_ylabel('Latitude ($\degree$)',
                        fontsize=FontPtSize,
                        fontweight='bold')
        csax.set_title(
            'Highest Concentration Reached Throughout\nAll Depths Anytime During Model Run',
            fontsize=FontPtSize + 2,
            fontweight='bold')
        cscbar = csfig.colorbar(cssurf)
        cscbar.set_label('Concentration (%s)' % units,
                         fontsize=FontPtSize,
                         fontweight='bold')
        cscbar.ax.tick_params(
            labelsize=FontPtSize)  #Do you want to add a tick size parameter?
        minX, maxX, minY, maxY = csax.axis()
        cbarfig = np.NaN
        #csfig.savefig(path4image, bbox_inches='tight')

    #TODO: Include in the function output a variable of the lat/lon coordinates of each contour and the concentration that the contour represents.
    #cssurf.levels
    #p=cssurf.collections[0].get_paths()[0] #First zero is nLevel, second zero is nContour for that level.
    #v = p.vertices
    #x = v[:,0]
    #y = v[:,1]
    #http://stackoverflow.com/questions/5666056/matplotlib-extracting-data-from-contour-lines

    #Calculate the total number of verticies used to define the contours.
    nVerticies = 0
    for i in np.arange(cssurf.collections.__len__()):
        nLevel = cssurf.collections[i].get_paths()
        for j in np.arange(nLevel.__len__()):
            nVerticies += len(nLevel[j])

    #TODO:  Create variables of the coordinates that form the contours of each level.
    #The following are some notes on how to properly get the coordinates of each contour (at each level).
    #i=0
    #nLevel = cssurf.collections[i].get_paths()
    #test=nLevel[6]
    #list(test.iter_segments())

    return (
        csfig,
        minX,
        maxX,
        minY,
        maxY,
        cbarfig,
    )


def make_pColorFig(Lon2D,
                   Lat2D,
                   masked_array,
                   units,
                   MinConc2Plot=None,
                   MaxConc2Plot=None,
                   FontPtSize=14,
                   isGoogleEarthKMZ=False):
    '''
    Create a figure with concentration plotted as a range of colors.
    Input:
        "Lon2D" is a two-dimensional rectangular grid of the longitudes
            corresponding to the concentration values in "masked_array".
        "Lat2D" is a two-dimensional rectangular grid of the latitudes
            corresponding to the concentration values in "masked_array".
        "masked_array" is a masked array of a rectangular grid of
            concentrations with zero values masked-out.
        "units" are the units of the concentration values in masked_array; either:
            "ppt" for parts per thousand (a.k.a. kilograms per cubic meter)
            "%" for percent
            "fraction" for fraction (decimal)
            "micro-g/l" for microgram per liter
            "mg/l" for milligrams per liter
            "ng/l" for nanograms per liter
            "ppb" for parts per billion
            "ppm" for parts per million
            "pptrillion" for parts per trillion
            "lb/ft^3" for pounds per cubic foot
            This "units" string is just used only to label the colorbar.
        "MinConc2Plot" is optional and is the lowest concentration that will be
            plotted with its own distinct color -- all values smaller than this
            threshold will be saturated with this color.  This is the lower
            limit of the color-bar.  Value must be in the same units as the
            concentrations in "masked_array".
        "MaxConc2Plot" is optional and is the highest concentration that will
            be plotted with its own distinct color -- all values larger than
            this threshold will be saturated with this color.  This is the
            upper limit of the color-bar.  Value must be in the same units as
            the concentrations in "masked_array".
        "FontPtSize" is optional and is the size of the font in points to use
            in labeling the plot.
        "isGoogleEarthKMZ" is optional and, if set to true, will return a
            figure without axes, ticks, or labels for use as an image overlay
            in Google Earth.
    Output:
        "pcfig" is an instance of a pseudocolor plot that can be used to generate
            or save the graphic using MatPlotLib.  See description of
            "isGoogleEarthKMZ" input.
        "minX" is the lower limit of the x (longitude) axis.
        "maxX" is the upper limit of the x (longitude) axis.
        "minY" is the lower limit of the y (latitude) axis.
        "maxY" is the upper limit of the y (latitude) axis.
        "cbarfig" is a separate figure instance that just contains the colorbar
            that goes with "pcfig".  A figure instance is returned if
            "isGoogleEarthKMZ" is True; otherwise a NaN is returned.
    Example:
        pcfig, minX, maxX, minY, maxY, cbarfig = plot_results.make_pColorFig( Lon2D, Lat2D, masked_array, MinConc2Plot, MaxConc2Plot, FontPtSize=14, isGoogleEarthKMZ=False )
    '''

    if isGoogleEarthKMZ:

        #Make a figure with without axes, ticks, or labels for use as an image overlay in Google Earth
        pcfig = plt.figure(frameon=False)
        pcax = pcfig.add_axes((0, 0, 1, 1))
        pcsurf = pcax.pcolor(Lon2D,
                             Lat2D,
                             masked_array,
                             shading='gouraud',
                             edgecolors='none')
        pcsurf.set_clim(MinConc2Plot, MaxConc2Plot)
        minX, maxX, minY, maxY = pcax.axis()
        pcax.set_axis_off()
        #pcfig.savefig(path4image, transparent=True, bbox_inches=0, pad_inches=0)

        #Create a seperate figure instance that just has the colorbar that goes with the pcolor plot.
        cbarfig = plt.figure()
        pccbar = cbarfig.colorbar(pcsurf)
        pccbar.set_label('Concentration (%s)' % units,
                         fontsize=FontPtSize,
                         fontweight='bold')
        pccbar.ax.tick_params(
            labelsize=FontPtSize)  #Do you want to add a tick size parameter?
        cbarfig.delaxes(cbarfig.gca())
        #cbarfig.savefig(path4cbar, transparent=False, bbox_inches='tight')

    else:

        #Make a complete figure, including an embedded colorbar (in which case the "cbarfig" output will likely be of no use).
        pcfig = plt.figure()
        pcax = pcfig.add_subplot(1, 1, 1)
        pcsurf = pcax.pcolor(Lon2D,
                             Lat2D,
                             masked_array,
                             shading='interp',
                             edgecolors='none')
        pcsurf.set_clim(MinConc2Plot,
                        MaxConc2Plot)  #Values outside range are saturated.
        pcax.ticklabel_format(style='plain', useOffset=False, axis='both')
        pcax.tick_params(labelsize=FontPtSize)
        pcax.set_xlabel('Longitude ($\degree$)',
                        fontsize=FontPtSize,
                        fontweight='bold')
        pcax.set_ylabel('Latitude ($\degree$)',
                        fontsize=FontPtSize,
                        fontweight='bold')
        pcax.set_title(
            'Highest Concentration Reached Throughout\nAll Depths Anytime During Model Run',
            fontsize=FontPtSize + 2,
            fontweight='bold')
        pccbar = pcfig.colorbar(pcsurf)
        pccbar.set_label('Concentration (%s)' % units,
                         fontsize=FontPtSize,
                         fontweight='bold')
        pccbar.ax.tick_params(
            labelsize=FontPtSize)  #Do you want to add a tick size parameter?
        minX, maxX, minY, maxY = pcax.axis()
        #pcfig.savefig(path4image, bbox_inches='tight')
        cbarfig = np.NaN

    return (
        pcfig,
        minX,
        maxX,
        minY,
        maxY,
        cbarfig,
    )


#==============================================================================
# Section #4:  File output functions.
#==============================================================================


def write_KMZ(path4KMZ, ColorBarFig, OverlayFig, WestLonLim, EastLonLim,
              SouthLatLim, NorthLatLim):
    '''
    Create a Google Earth KMZ file with an overlay of the associated color-bar.
    Input:
        "path4KMZ" is the full path of where the KMZ file is to be saved.
        "ColorBarFig" is the color-bar figure object, that provides a legend of
            the concentrations plotted in the overlay figure.
        "OverlayFig" is the main figure object, which shows a 2-D spatial
            distribution of concentration.
        "WestLonLim" is the western longitude extent of the overlay figure, in
            units of decimal degrees.
        "EastLonLim" is the eastern longitude extent of the overlay figure, in
            units of decimal degrees.
        "SouthLatLim" is the southern latitude extent of the overlay figure, in
            units of decimal degrees.
        "NorthLatLim" is the northern latitude extent of the overlay figure, in
            units of decimal degrees.
    Output:
        A KMZ file saved at the location specified.
    Example:
        plot_results.write_KMZ( 'C:\\Users\brian.zelenke\Documents\GNOME\Tests\GridPlume\ContinuousReleaseUniformCurrent\20130328_Zelenke_GridPlumeOutput_UniformSpillNoDiff.kmz', cbfig, fig, -92.05, -91.7, 26.998, 27.002 )
    '''

    #Specify subdirectory where image files will be written.
    filesdir = path4KMZ[:]
    filesdir += "/files"

    #Path for image overlay (a.k.a. tile); in the KML "files" directory.
    path4fig = filesdir[:]
    path4fig += "/overlay.png"

    #Path for colorbar legend; in the KML "files" directory.
    path4cbar = filesdir[:]
    path4cbar += "/ColorBar.png"

    #Path for the KML document.
    path4KML = path4KMZ[:]
    path4KML += "/doc.kml"

    #Path for temporary .zip file.
    tempzip = path4KMZ[:]
    tempzip += ".zip"
    #Delete it if it exists... only likely if this plotting routine got interrupted.
    if os.path.exists(tempzip):
        try:
            shutil.rmtree(tempzip)
        except:
            os.remove(tempzip)

    #Create a directory whose name includes the .kmz extension.
    if os.path.exists(path4KMZ):
        try:
            shutil.rmtree(path4KMZ)
        except:
            #Different exceptions are raised by shutil depending on the operating system, so don't make this exception-specific.
            os.remove(path4KMZ)
    os.makedirs(filesdir)

    #Write the image overlay (a.k.a. tile) and colorbar (a.k.a legend) into the KML "files" directory.
    OverlayFig.savefig(path4fig, transparent=True, bbox_inches=0, pad_inches=0)
    ColorBarFig.savefig(path4cbar, transparent=False, bbox_inches='tight')

    #Create ISO 8601 formatted string of the UTC time (March 1, 2000 15:45:17 UTC = 20000301T154517 per format string 'yyyymmddTHHMMSS').  If you wanted to include microseconds for more precision, remove the "[0:17]".
    now = datetime.datetime.utcnow().isoformat().replace(":", "")[0:17]

    #Write the KML code describing how to plot the tile and colorbar images into a file.
    fKML = open(path4KML, 'w')
    fKML.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    fKML.write(
        '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n'
    )

    fKML.write('<Folder>\n')
    fKML.write('\t<name>GridPlume</name>\n')
    fKML.write(
        '\t<description>Highest Concentration Reached Throughout All Depths Anytime During Model Run.\n'
    )
    fKML.write('Created by GridPlume at %s.</description>\n' %
               now)  #Make the description of the overlay include a time-stamp.

    fKML.write('\t<ScreenOverlay>\n')
    fKML.write('\t\t<name>Color-bar</name>\n')
    fKML.write('\t\t<Icon>\n')
    fKML.write('\t\t\t<href>files/ColorBar.png</href>\n')
    fKML.write('\t\t</Icon>\n')
    fKML.write(
        '\t\t<overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n')
    fKML.write(
        '\t\t<screenXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n')
    fKML.write(
        '\t\t<rotationXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>\n'
    )
    fKML.write('\t\t<size x="-1" y="-1" xunits="pixels" yunits="pixels"/>\n')
    fKML.write('\t</ScreenOverlay>\n')

    fKML.write('\t<GroundOverlay>\n')
    fKML.write('\t\t<name>MaxConcentration</name>\n')
    fKML.write('\t\t<Icon>\n')
    fKML.write('\t\t\t<href>files/overlay.png</href>\n')
    fKML.write('\t\t</Icon>\n')
    fKML.write('\t\t<LatLonBox>\n')
    fKML.write('\t\t\t<north>%g</north>\n' % NorthLatLim)
    fKML.write('\t\t\t<south>%g</south>\n' % SouthLatLim)
    fKML.write('\t\t\t<east>%g</east>\n' % EastLonLim)
    fKML.write('\t\t\t<west>%g</west>\n' % WestLonLim)
    fKML.write('\t\t</LatLonBox>\n')
    fKML.write('\t</GroundOverlay>\n')

    fKML.write('</Folder>\n')
    fKML.write('</kml>')
    fKML.close()

    #Zip the KML directory, thus creating a KMZ.
    shutil.make_archive(path4KMZ, format="zip", root_dir=path4KMZ)

    #Delete the KML directory, leaving only the KMZ.
    shutil.rmtree(path4KMZ)

    #Take off the ".zip" file extension, leaving just the ".kmz".
    zipdir = path4KMZ[:]
    zipdir += ".zip"
    os.rename(zipdir, path4KMZ)


#==============================================================================
# Section #5:  Main plotting functions.
#==============================================================================


def plot_max_conc_everywhere(path2file,
                             units,
                             path4image=None,
                             MinConc2Plot=None,
                             MaxConc2Plot=None,
                             Just3Colors=False,
                             LowConcUpperLim=None,
                             HighConcLowerLim=None,
                             FontPtSize=14,
                             levels=None):
    #TODO: Include a "path4file" variable to invoke the NetCDF/CSV file output options of the above "max_conc_grids" function.
    '''
    Plot the maximum concentration modeled at each grid-point during the
    sampling time, irrespective of depth.  Optional output can be an image file
    or a Google Earth KMZ file.
    Input:
        "path2file" is the full path to a NetCDF file output by grid_plume.py.
            The NetCDF file must contain the following variables in the units given:
                "lat" in units of decimal degrees north;
                "lon" in units of decimal degrees east;
                "depth" in units of meters;
                "time" in units of seconds since a date-time (e.g., Unix time)
                "concentration" in units of kilograms per cubic meter.
         "units" are the units concentration will be plotted in; either:
            "ppt" for parts per thousand (a.k.a. kilograms per cubic meter)
            "%" for percent
            "fraction" for fraction (decimal)
            "micro-g/l" for microgram per liter
            "mg/l" for milligrams per liter
            "ng/l" for nanograms per liter
            "ppb" for parts per billion
            "ppm" for parts per million
            "pptrillion" for parts per trillion
            "lb/ft^3" for pounds per cubic foot
        "path4image" is optional and is a string containing a path to the
            file-name of the image you want to save. The file-name can end in
            either:
                -- any common image extension (".png" is recommended, but
                ".jpg" and ".gif" are available too); or
                -- ".kmz", in which case a Google Earth KMZ file of the image
                    and with its color-bar will be saved.
        "MinConc2Plot" is optional and is the lowest concentration that will be
            plotted with its own distinct color -- all values smaller than this
            threshold will be saturated with this color.  This is the lower
            limit of the color-bar.
        "MaxConc2Plot" is optional and is the highest concentration that will
            be plotted with its own distinct color -- all values larger than
            this threshold will be saturated with this color.  This is the
            upper limit of the color-bar.
        "Just3Colors" is optional and either a logical True or False (default
            is false) where True will plot the concentration color-scale using
            only three contours of "high-medium-low" colors and False will plot
            the concentration color-scale using a rainbow of continuous colors.
        "LowConcUpperLim" is optional and will be used only if "Just3Colors" is
            True.  LowConcUpperLim specifies the upper limit of the values to
            consider "low" in concentration.  All values between MinConc2Plot
            up through LowConcUpperLim will be shaded blue.
        "HighConcLowerLim" is optional and will be used only if "Just3Colors"
            is true.  All values between HighConcLowerLim up through
            MaxConc2Plot will be shaded green.
        "TransparentBackground" is optional and, if True, will make the
            background of the figure transparent (handy for overlaying the
            image on maps, etc.).
        "FontPtSize" is optional and is the size of the font in points to use
            in labeling the plot.  The title will be two point sizes larger.
    Output:
        A figure shown on-screen and, if specified, an image or Google Earth
            file.
        "fig" is the main figure object, which shows a 2-D spatial distribution
            of concentration.
        "cbfig" is the color-bar figure object, that provides a legend of the
            concentrations plotted in the overlay figure.

    Examples:
        plot_results.plot_max_conc_everywhere( 'C:/Users/brian.zelenke/Documents/GNOME/UXO/GlenWatabayashi/UXO_ROMS/Zelenke_test/20120315_GridPlumeOutput3Days10kg1000les_uh.nc', 'ppb', MinConc2Plot=0, MaxConc2Plot=1 )
        plot_results.plot_max_conc_everywhere( 'C:/Users/brian.zelenke/Documents/GNOME/Tests/GNOME_3d_test/2013-03-14_PresentationToDebbiePayton/20130312_GridPlumeOutput_Zelenke_All_spills_500m.nc', 'ppm', MinConc2Plot=0, MaxConc2Plot=50 )
    '''

    Lon2D, Lat2D, MaxConcGrid = max_conc_grids(path2file, units, levels=levels)

    #Mask-out areas with no concentration values.
    masked_array = np.ma.array(MaxConcGrid, mask=np.isnan(MaxConcGrid))

    #Calculate color limit values the user didn't specify.
    if MinConc2Plot == None:
        MinConc2Plot = np.nanmin(MaxConcGrid)
    if MaxConc2Plot == None:
        MaxConc2Plot = np.nanmax(MaxConcGrid)
    MinConc2Plot, MaxConc2Plot = np.sort([MinConc2Plot,
                                          MaxConc2Plot])  #Quality control.

    #Determine if output requested is a Google Earth KMZ file.
    if path4image == None:
        isGoogleEarthKMZ = False
    elif path4image[-4:].lower() == ".kml":
        #If user accidently specifies a Google Earth KML file, correct it to KMZ.
        path4image = path4image[0:-4]
        path4image += '.kmz'
        isGoogleEarthKMZ = True
    else:
        isGoogleEarthKMZ = path4image[-4:].lower() == ".kmz"

    #Create figure instances.
    if Just3Colors:

        #(Lightly) quality control any user-specified color axis ticks and
        #calculate any missing values.
        levels = np.linspace(np.nanmin(MinConc2Plot), np.nanmax(MaxConc2Plot),
                             4)
        if LowConcUpperLim == None:
            LowConcUpperLim = levels[1]
        if HighConcLowerLim == None:
            HighConcLowerLim = levels[2]
        LowConcUpperLim, HighConcLowerLim = np.sort(
            [LowConcUpperLim, HighConcLowerLim])
        if LowConcUpperLim > MinConc2Plot and LowConcUpperLim < MaxConc2Plot:
            levels[1] = LowConcUpperLim
        if HighConcLowerLim > MinConc2Plot and HighConcLowerLim < MaxConc2Plot:
            levels[2] = HighConcLowerLim

        fig, minX, maxX, minY, maxY, cbfig = make_3ColorFig(
            Lon2D, Lat2D, masked_array, units, levels, FontPtSize,
            isGoogleEarthKMZ)

    else:
        fig, minX, maxX, minY, maxY, cbfig = make_pColorFig(
            Lon2D, Lat2D, masked_array, units, MinConc2Plot, MaxConc2Plot,
            FontPtSize, isGoogleEarthKMZ)

    #Write-out the figure as either an image file or Google Earth KMZ.
    if isGoogleEarthKMZ:
        write_KMZ(path4image, cbfig, fig, minX, maxX, minY, maxY)
    elif path4image != None:
        fig.savefig(path4image, bbox_inches='tight')

    plt.show(
        fig
    )  #Show the figure on the screen.  If a KMZ was specified, this will just be the tile sans labels.

    return (
        fig,
        cbfig,
    )


def plot_conc_at_point(path2file,
                       units,
                       PointLat,
                       PointLon,
                       path4image=None,
                       FontPtSize=14):
    #TODO: Include a "path4file" variable to invoke the NetCDF/CSV file output options of the above "point_conc" function.
    '''
    Plot the time-series of the highest concentration modeled at any depth at
    the grid location closest to the coordinate given, and save to an image
    file if requested.
    Input:
        "path2file" is the full path to a NetCDF file output by grid_plume.py.
            The NetCDF file must contain the following variables in the units given:
                "lat" in units of decimal degrees north;
                "lon" in units of decimal degrees east;
                "depth" in units of meters;
                "time" in units of seconds since a date-time (e.g., Unix time)
                "concentration" in units of kilograms per cubic meter.
        "units" are the units concentration will be plotted in; either:
            "ppt" for parts per thousand (a.k.a. kilograms per cubic meter)
            "%" for percent
            "fraction" for fraction (decimal)
            "micro-g/l" for microgram per liter
            "mg/l" for milligrams per liter
            "ng/l" for nanograms per liter
            "ppb" for parts per billion
            "ppm" for parts per million
            "pptrillion" for parts per trillion
            "lb/ft^3" for pounds per cubic foot
        "PointLat" is the latitude coordinate in decimal degrees (+=northern
            hemisphere, -=southern hemisphere) of your location of interest
            that falls within the area modeled with GridPlume.
        "PointLon" is the longitude coordinate in decimal degrees (+=eastern
            hemisphere, -=western hemisphere) of your location of interest that
            falls within the area modeled with GridPlume.
        "path4image" is optional and is a string containing a path to the
            file-name of the image you want to save. The file-name can end in
            any common image extension (".png" is recommended, but ".jpg"
            and ".gif" are available too).
        "FontPtSize" is optional and is the size of the font in points to use
            in labeling the plot.  The title will be two point sizes larger.
    Output:
        A figure shown on-screen and, if specified, an image file.
        "fig" is the figure object of the image shown on-screen.
    Example:
        plot_results.plot_conc_at_point( 'C:/Users/brian.zelenke/Documents/GNOME/UXO/GlenWatabayashi/UXO_ROMS/Zelenke_test/20120315_GridPlumeOutput3Days10kg1000les_uh.nc', 'ppb', 21.453471, -158.204816 )
    '''

    #Get time-series of maximum concentration (irrespective of depth) at
    #grid-point nearest coordinate given.
    MaxConcAtPt, time, pLat, pLon, _ = point_conc(path2file, units, PointLat,
                                                  PointLon)

    tHours = time / (
        60 * 60
    )  #Convert time from seconds elapsed (Unix time) to hours elapsed.

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    l = ax.plot(tHours - tHours[0], MaxConcAtPt, 'k-', linewidth=2)
    ax.set_xlim(0, np.max(tHours - tHours[0]))
    ax.tick_params(labelsize=FontPtSize)
    ax.set_xlabel('Hours from Release Time',
                  fontsize=FontPtSize,
                  fontweight='bold')
    ax.set_ylabel('Concentration (%s)' % units,
                  fontsize=FontPtSize,
                  fontweight='bold')
    ax.set_title(
        r'Time-series of Highest Concentration within Water Column at %g$^{\circ}$, %g$^{\circ}$'
        % (pLat, pLon),
        fontsize=FontPtSize + 2,
        fontweight='bold')

    #Write-out the figure as either an image file
    if path4image != None:
        fig.savefig(path4image, bbox_inches='tight')

    plt.show(fig)  #Show the figure on the screen.

    return (fig)


def plot_max_conc_in_plume(path2file, units, path4image=None, FontPtSize=14):
    #TODO: Include a "path4file" variable to invoke the NetCDF/CSV file output options of the above "max_conc_in_plume" function.
    '''
    Plot the time-series of the highest concentration modeled at any depth
    anywhere within the given grid, and save to an image file if requested.
    Input:
        "path2file" is the full path to a NetCDF file output by grid_plume.py.
            The NetCDF file must contain the following variables in the units given:
                "lat" in units of decimal degrees north;
                "lon" in units of decimal degrees east;
                "depth" in units of meters;
                "time" in units of seconds since a date-time (e.g., Unix time)
                "concentration" in units of kilograms per cubic meter.
        "units" are the units concentration will be plotted in; either:
            "ppt" for parts per thousand (a.k.a. kilograms per cubic meter)
            "%" for percent
            "fraction" for fraction (decimal)
            "micro-g/l" for microgram per liter
            "mg/l" for milligrams per liter
            "ng/l" for nanograms per liter
            "ppb" for parts per billion
            "ppm" for parts per million
            "pptrillion" for parts per trillion
            "lb/ft^3" for pounds per cubic foot
        "path4image" is optional and is a string containing a path to the
            file-name of the image you want to save. The file-name can end in
            any common image extension (".png" is recommended, but ".jpg"
            and ".gif" are available too).
        "FontPtSize" is optional and is the size of the font in points to use
            in labeling the plot.  The title will be two point sizes larger.
    Output:
        A figure shown on-screen and, if specified, an image file.
        "fig" is the figure object of the image shown on-screen.
    Example:
        plot_results.plot_max_conc_in_plume( 'C:/Users/brian.zelenke/Documents/GNOME/UXO/GlenWatabayashi/UXO_ROMS/Zelenke_test/20120315_GridPlumeOutput3Days10kg1000les_uh.nc', 'ppb' )
    '''

    #Get the time-series of maximum concentration (irrespective of depth) within the plume(s).
    MaxConcSeries, time = max_conc_in_plume(path2file, units, levels)

    tHours = time / (
        60 * 60
    )  #Convert time from seconds elapsed (Unix time) to hours elapsed.

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    l = ax.plot(tHours - tHours[0], MaxConcSeries, 'k-', linewidth=2)
    ax.set_xlim(0, np.max(tHours - tHours[0]))
    ax.tick_params(labelsize=FontPtSize)
    ax.set_xlabel('Hours from Release Time',
                  fontsize=FontPtSize,
                  fontweight='bold')
    ax.set_ylabel('Concentration (%s)' % units,
                  fontsize=FontPtSize,
                  fontweight='bold')
    ax.set_title(
        r'Time-series of Highest Concentration within Water Column Throughout Model Domain',
        fontsize=FontPtSize + 2,
        fontweight='bold')

    #Write-out the figure as either an image file
    if path4image != None:
        fig.savefig(path4image, bbox_inches='tight')

    plt.show(fig)  #Show the figure on the screen.

    return (fig)
