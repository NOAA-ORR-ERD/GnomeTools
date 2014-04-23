from __future__ import division #Change the / operator to ensure true division throughout (Zelenke).
#import operator #Commented-out import of unused module (Zelenke).
import pyproj
import numpy as np
import netCDF4 as nc


def init_save_grid(
    dataset,
    grid_shape,
    times,
    time_units,
    projection,
    x_min,
    x_max,
    y_min,
    y_max,
    z_min,
    z_max,
    units,
):
    """
    Prepare a netCDF4 dataset for saving a concentration grid. This
    preparation involves creating the appropriate dimensions, variables, and
    attributes.

    Return the prepared grid variable object so that the caller can write grid
    data to it. The caller is also responsible for closing the dataset when
    done.

    :param dataset: open writeable netCDF4 dataset
    :type dataset: netCDF.Dataset
    :param grid_shape: shape of concentration grid, e.g. grid.shape
    :type grid_shape: tuple
    :param times: a time for each time step
    :type times: array of datetime objects
    :param time_units: netCDF4 units to use for the times variable
    :type time_units: str
    :param projection: projection used for the following coordinates
    :type projection: pyproj.Proj
    :param x_min: the minimum x coordinate of the grid in meters
    :type x_min: float
    :param x_max: the maximum x coordinate of the grid in meters
    :type x_max: float
    :param y_min: the minimum y coordinate of the grid in meters
    :type y_min: float
    :param y_max: the maximum y coordinate of the grid in meters
    :type y_max: float
    :param z_min: the minimum z coordinate of the grid in meters
    :type z_min: float
    :param z_max: the maximum z coordinate of the grid in meters
    :type z_max: float
    :param units: units to use for concentration grid
    :type units: str
    """
    # Convert meters to latitude and longitude.
    lat_long = pyproj.Proj( "+proj=latlong" )

    x = np.linspace( x_min, x_max, grid_shape[ 2 ] )
    y = np.linspace( y_min, y_max, grid_shape[ 1 ] )
    z = np.linspace( z_min, z_max, grid_shape[ 0 ] )

    ( longitude, _ ) = pyproj.transform( projection, lat_long, x, x )
    ( _, latitude ) = pyproj.transform( projection, lat_long, y, y )
    

    dataset.createDimension( "depth", grid_shape[ 0 ] )
    dataset.createDimension( "lon", grid_shape[ 2 ] )
    dataset.createDimension( "lat", grid_shape[ 1 ] )
    dataset.createDimension( "time", times.shape[ 0 ] )

    nc_longitude = dataset.createVariable(
        "lon",
        np.float32,
        "lon",
        zlib = True,
    )
    nc_longitude.long_name = "Longitude"
    nc_longitude.units = "degrees_east"
    nc_longitude.standard_name = "longitude"
    nc_longitude[ : ] = longitude

    nc_latitude = dataset.createVariable(
        "lat",
        np.float32,
        "lat",
        zlib = True,
    )
    nc_latitude.long_name = "Latitude"
    nc_latitude.units = "degrees_north"
    nc_latitude.standard_name = "latitude"
    nc_latitude[ : ] = latitude

    nc_depth = dataset.createVariable(
        "depth",
        np.float32,
        "depth",
        zlib = True,
    )
    nc_depth.long_name = "Depth"
    nc_depth.units = "meters"
    nc_depth.positive = "up"
    nc_depth.standard_name = "depth"
    nc_depth[ : ] = z

    nc_time = dataset.createVariable(
        "time",
        np.float32,
        "time",
        zlib = True,
    )
    nc_time.long_name = "Time"
    nc_time.units = time_units
    nc_time.standard_name = "time"
    nc_time[ : ] = nc.date2num(
        times,
        units = time_units,
    )
    #The lat, lon, depth, and mass coming out of GNOME are all single precision.
    #The concentration calculation doesn't add any significant figures beyond
    #that, so no improvement seems gained for the extra memory.  Thus the
    #double-precision fill value has been commented-out below.
    nc_grid = dataset.createVariable(
        "concentration",
        np.float32,
        ( "time", "depth", "lat", "lon" ),
        zlib = True,
        #fill_value = nc.default_fillvals[ "f4" ], #(Zelenke)
    )
    #Back when grid_plume.py created a NetCDF dataset with the format
    #NETCDF3_CLASSIC, until this point the variables didn't actually seem to be
    #written-out to the NetCDF file.  Once the character strings below were
    #specified, however, the NetCDF file (as seen by the OS) ballooned to
    #full-size.  Not sure why specifying these strings triggers the
    #writing/ballooning of the NetCDF file.  (Zelenke)
    nc_grid.long_name = "Mass concentration of total hydrocarbons in seawater"
    nc_grid.units = units
    nc_grid.standard_name = "mass_concentration_of_total_hydrocarbons_in_sea_water"

    setattr( dataset, "Conventions", "CF-1.4" )
    setattr( dataset, "title", "Oil Plume Concentration Grid" )
    setattr( dataset, "institution", "NOAA/NOS/OR&R/ERD/TSSB" )
    setattr( dataset, "references", "brian.zelenke@noaa.gov" )

    return nc_grid


def init_flat_save_grid(
    dataset,
    grid_shape,
    times,
    time_units,
    projection,
    x_min,
    x_max,
    y_min,
    y_max,
    units,
):
    """
    Prepare a netCDF4 dataset for saving a concentration grid. This
    preparation involves creating the appropriate dimensions, variables, and
    attributes.

    Return the prepared grid variable object so that the caller can write grid
    data to it. The caller is also responsible for closing the dataset when
    done.

    :param dataset: open writeable netCDF4 dataset
    :type dataset: netCDF.Dataset
    :param grid_shape: shape of flattened concentration grid, e.g. grid.shape
    :type grid_shape: tuple
    :param times: a time for each time step
    :type times: array of datetime objects
    :param time_units: netCDF4 units to use for the times variable
    :type time_units: str
    :param projection: projection used for the following coordinates
    :type projection: pyproj.Proj
    :param x_min: the minimum x coordinate of the grid in meters
    :type x_min: float
    :param x_max: the maximum x coordinate of the grid in meters
    :type x_max: float
    :param y_min: the minimum y coordinate of the grid in meters
    :type y_min: float
    :param y_max: the maximum y coordinate of the grid in meters
    :type y_max: float
    :param units: units to use for concentration grid
    :type units: str
    """
    # Convert meters to latitude and longitude.
    lat_long = pyproj.Proj( "+proj=latlong" )

    x = np.linspace( x_min, x_max, grid_shape[ 2 ] )
    y = np.linspace( y_min, y_max, grid_shape[ 1 ] )

    ( longitude, _ ) = pyproj.transform( projection, lat_long, x, x )
    ( _, latitude ) = pyproj.transform( projection, lat_long, y, y )

    dataset.createDimension( "lon", grid_shape[ 2 ] )
    dataset.createDimension( "lat", grid_shape[ 1 ] )
    dataset.createDimension( "time", times.shape[ 0 ] )

    nc_longitude = dataset.createVariable(
        "lon",
        np.float32,
        "lon",
        zlib = True,
    )
    nc_longitude.long_name = "Longitude"
    nc_longitude.units = "degrees_east"
    nc_longitude.standard_name = "longitude"
    nc_longitude[ : ] = longitude

    nc_latitude = dataset.createVariable(
        "lat",
        np.float32,
        "lat",
        zlib = True,
    )
    nc_latitude.long_name = "Latitude"
    nc_latitude.units = "degrees_north"
    nc_latitude.standard_name = "latitude"
    nc_latitude[ : ] = latitude

    nc_time = dataset.createVariable(
        "time",
        np.float32,
        "time",
        zlib = True,
    )
    nc_time.long_name = "Time"
    nc_time.units = time_units
    nc_time.standard_name = "time"
    nc_time[ : ] = nc.date2num(
        times,
        units = time_units,
    )

    nc_grid = dataset.createVariable(
        "concentration",
        np.float32,
        ( "time", "lat", "lon" ),
        zlib = True,
        fill_value = nc._default_fillvals[ "f4" ],
    )
    nc_grid.long_name = "Mass concentration of total hydrocarbons in seawater"
    nc_grid.units = units
    nc_grid.standard_name = "mass_concentration_of_total_hydrocarbons_in_sea_water"

    setattr( dataset, "Conventions", "CF-1.4" )
    setattr( dataset, "title", "Oil Plume Concentration Grid" )
    setattr( dataset, "institution", "NOAA/NOS/ERD/TSSB" )
    setattr( dataset, "references", "brian.zelenke@noaa.gov" )

    return nc_grid


def flatten_concentration(
    file_grid,
    flat_file_grid,
    grid_z_min,
    grid_z_max,
):
    bucket_z_delta = ( grid_z_max - grid_z_min ) / file_grid.shape[ 1 ]

    for time_index in range( file_grid.shape[ 0 ] ):
        flat_file_grid[ time_index ] = \
            file_grid[ time_index ].sum( axis = 0 ) / bucket_z_delta
