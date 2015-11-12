#!/usr/bin/env python

from __future__ import division #Change the / operator to ensure true division throughout (Zelenke).
print "Running grid_plume."
print "Running with GNOME particle data extension"

import os
#import sys #Commented-out import of unused module (Zelenke).
import time
from optparse import OptionParser
#import itertools #Commented-out import of unused module (Zelenke).
import numpy as np

import gridplume
from gridplume import filter
from gridplume import create_grid

import netCDF4 as nc
import pyproj
#import BuildSample 
#import nc_particles
#import TestRead

DEFAULT_GRID_STEPS = 64
DEFAULT_SIGMA_0 = 200.0 # meters
DEFAULT_D_X = 0.01 # m^2/s
DEFAULT_D_Z = 0.001 # m^2/s

DEFAULT_ALPHA = 0.5
DEFAULT_MASS = 338.56
DEFAULT_WATER_DENSITY = 1000.0
DEFAULT_OUTPUT_UNITS = "ppm"

OUTPUT_UNITS_MAP = dict(
    density = "kg/m^3",
    flat = "kg/m^2",
    ratio = "kg/kg",
    ppm = "ppm",
    ppb = "ppb",
)

# See http://www.unidata.ucar.edu/software/netcdf/docs/netcdf/NetCDF-Classic-Format-Limitations.html#NetCDF-Classic-Format-Limitations (doesn't much matter now, as GridPlume now outputs a NetCDF4 file.)
#MAX_NETCDF_CLASSIC_RECORD_SIZE = 2 ** 32 - 4 #Old NetCDF3 Classic format limitation.
MAX_NETCDF_CLASSIC_RECORD_SIZE = np.inf #No maximum now that we're using NetCDF4.
MIN_BUCKET_COUNT_X = 32
MIN_BUCKET_COUNT_Y = 32
MIN_BUCKET_COUNT_Z = 8
BOGUS_PARTICLE_FLAG = -99999.0

USAGE = "usage: %prog [options] particlefile [particlefile, ...] outputfile"
MINIMUM_ARGS_COUNT = 2


def parse_args():
    parser = OptionParser( usage = USAGE )

    parser.add_option(
        "-s", "--sigma_0", dest = "sigma_0",
        help = "initial x dimension of the parcels in meters (%s)" %
            DEFAULT_SIGMA_0,
        default = DEFAULT_SIGMA_0,
        type = "float",
    )
    parser.add_option(
        "-x", "--diff_x", dest = "D_x",
        help = " horizontal diffusion coeff in m^2/s (%s)" %
            DEFAULT_D_X,
        default = DEFAULT_D_X,
        type = "float",
    )
    parser.add_option(
        "-z", "--diff_z", dest = "D_z",
        help = " vertical diffusion coeff in m^2/s (%s)" %
            DEFAULT_D_Z,
        default = DEFAULT_D_Z,
        type = "float",
    )
    parser.add_option(
        "-a", "--alpha", dest = "alpha",
        help = "exponent of spreading (%s)" % DEFAULT_ALPHA,
        default = DEFAULT_ALPHA,
        type = "float",
    )
    parser.add_option(
        "-m", "--mass_file", dest = "mass_file",
        help = "file to calculate particle mass from (mass defaults to %s kg if no file is given)" % DEFAULT_MASS,
        type = "string",
    )
    parser.add_option(
        "-w", "--water-density", dest = "water_density",
        help = 'density of the water in kg/m^3, ignored if output units is "density" (%s)' % DEFAULT_WATER_DENSITY,
        default = DEFAULT_WATER_DENSITY,
        type = "float",
    )
    parser.add_option(
        "-u", "--output_units", dest = "output_units",
        help = 'units to use in output concentration grid: one of "density" (meaning kg/m^3), "flat" (meaning one file with kg/m^3 and a separate file with kg/m^2 and depth vertically integrated), "ratio" (meaning kg oil per kg water), "ppm" (meaning parts per million), or "ppb" (meaning parts per billion) (defaults to "%s")' % DEFAULT_OUTPUT_UNITS,
        default = DEFAULT_OUTPUT_UNITS,
        type = "choice",
        choices = OUTPUT_UNITS_MAP.keys(),
    )
    parser.add_option(
        "-g", "--grid_file", dest = "grid_file",
        help = 'file with grid specification: it should specify: lon_min, lon_max, lat_min, lat_max, z_min, z_max, num_lon, num_lat, num_z (in lat-long) in python format if not specified, it will be computed to fit particles',
        default = None,
        type = "string",
    )

    ( options, args ) = parser.parse_args()

    if len( args ) < MINIMUM_ARGS_COUNT:
        parser.print_help()
        raise TypeError("grid_plume takes at least %g arguments (%g given)" % (MINIMUM_ARGS_COUNT, len(args)) )
        #sys.exit( 1 )

    return ( options, args )

def load_grid(filename, proj_lat_long, proj_meters):
    print "loading:", filename
    grid_params = {}
    execfile(filename, grid_params)
    lon_min = grid_params['lon_min']
    lon_max = grid_params['lon_max']
    lat_min = grid_params['lat_min']
    lat_max = grid_params['lat_max']

    # project from lat-long
    ( x_min, y_min ) = pyproj.transform( proj_lat_long, proj_meters, lon_min, lat_min )
    ( x_max, y_max ) = pyproj.transform( proj_lat_long, proj_meters, lon_max, lat_max )

    grid = np.zeros( (grid_params['num_z'],
                      grid_params['num_lat'],
                      grid_params['num_lon'])
                     , np.double)
    return(grid,
           x_min,
           x_max,
           y_min,
           y_max,
           grid_params['z_min'],
           grid_params['z_max'],
           )

def main():
    ( options, args ) = parse_args()
    particle_filename = args[ -2 ]
    output_file = args[ -1 ]

    print "Loading particles from particle file."

    if not os.path.isfile( particle_filename ):
        raise IOError( "The particle file \"%s\" does not exist.\nAborting." % particle_filename )
        #sys.exit( 1 )

    x, y, z, age, times, time_units, mass, mass_units = gridplume.load_particles(nc.Dataset( particle_filename, "r" )) #Added age and units of mass to output variables (Zelenke).

#==============================================================================
#     # check for all particles having the same mass
#     if np.allclose( mass[mass==mass], np.unique(mass[mass==mass])[0] ):
#         #Insist that all mass values be EXACTLY the same.
#         mass = np.array( np.mean(mass[mass==mass]), dtype=mass.dtype )
#     else:
#         raise ValueError( "The masses of the particles in file \"%s\" differ.\ngird_plume currently only supports particles of all the same mass.\nAborting." % particle_filename )
#==============================================================================


    #Added "if" statement below to ensure units of mass are kilograms (Zelenke).
    mass_units=mass_units.lower()
    if mass_units=='grams' or mass_units=='g' or mass_units=='gm':
        mass=mass/1000.0 #Convert g to kg.

    # Look for flagged bogus particles, and mark them as inf.
    # x[ np.equal( x, BOGUS_PARTICLE_FLAG ) ] = np.inf
    # y[ np.equal( y, BOGUS_PARTICLE_FLAG ) ] = np.inf
    # z[ np.equal( z, BOGUS_PARTICLE_FLAG ) ] = np.inf

    # Convert latitude and longitude values to meters.
    min_lon = np.nanmin( x )
    min_lat = np.nanmin( y )

    lat_long = pyproj.Proj( "+proj=latlong" )
    meters = pyproj.Proj(
        "+proj=merc +lon_0=%s +lat_ts=%s +x_0=0 +y_0=0 +units=m" % \
            ( min_lon, min_lat )
    )
    ( x, y ) = pyproj.transform( lat_long, meters, x, y )

    # convert to doubles:
    x = x.astype(np.float64)
    y = y.astype(np.float64)
    z = z.astype(np.float64)
    mass = mass.astype(np.float64)  #(Zelenke)
    age = age.astype(np.float64)   #(Zelenke)


    # Now that pyproj is through with the data, convert infs to NaNs.
    x[ np.isinf( x ) ] = np.nan
    y[ np.isinf( y ) ] = np.nan
    z[ np.isinf( z ) ] = np.nan
    mass[ np.isinf( mass ) ] = np.nan
    age[ np.isinf( age ) ] = np.nan

    # Make an ordered list of all the time steps from all the different
    # particle files, filtering out duplicates.
    seconds = gridplume.convert_times_to_seconds_since_start( times )

    #print "Initial particle mass: %s kg" % mass
    print "Initial masses of particles: %s kg" % np.unique( mass[ ~np.isnan(mass) ] ) #(Zelenke)
    #print "Initial total mass: %s kg" % ( len( x[ 0 ] ) * mass )
    print "Initial total mass: %s kg" % np.nansum( mass[0,:] ) #(Zelenke)

    ## compute the Coefficents:
    D_x = options.D_x
    D_z = options.D_z
    alpha = options.alpha
    C_x = 2*D_x
    C_y = C_x
    t_0 = options.sigma_0**(1/alpha) / (2*D_x)

    print "t_0: %f sec, %f hours"%( t_0,  t_0/3600 )

    C_z = 2*D_z

    if options.grid_file is None:
        print "Computing concentration grid."
        (
            grid,
            grid_x_min,
            grid_x_max,
            grid_y_min,
            grid_y_max,
            grid_z_min,
            grid_z_max,
        ) = create_grid.create_grid(
            x,
            y,
            z,
            seconds,
            C_x,
            C_y,
            C_z,
            t_0,
            options.alpha,
            MIN_BUCKET_COUNT_X,
            MIN_BUCKET_COUNT_Y,
            MIN_BUCKET_COUNT_Z,
            MAX_NETCDF_CLASSIC_RECORD_SIZE,
            # We're writing 32-bit floats out to file, so bucket size is 4 bytes.
            bucket_size = 4,
        )
    else:
        print "Loading grid definition from file: %s"%options.grid_file
        (grid,
         grid_x_min,
         grid_x_max,
         grid_y_min,
         grid_y_max,
         grid_z_min,
         grid_z_max,
        ) = load_grid(options.grid_file, lat_long, meters)

    print "Grid shape: %s (ZxYxX)" % "x".join( [ str( s ) for s in grid.shape ] )
    print "zmax:", grid_z_max, "zmin:", grid_z_min
    print "delta_x = ", (grid_x_max - grid_x_min)/ ( grid.shape[2] - 1.0 )
    print "delta_y = ", (grid_y_max - grid_y_min)/ ( grid.shape[1] - 1.0 )
    print "delta_z = ", (grid_z_max - grid_z_min)/ ( grid.shape[0] - 1.0 )

    # Write grid out to a new netCDF file.
    if options.output_units == "flat":
        print "Creating density concentration grid file for writing: %s." % \
            output_file
    else:
        print "Creating %s concentration grid file for writing: %s." % \
            ( options.output_units , output_file )
    #dataset = nc.Dataset( output_file, "w", format = "NETCDF3_CLASSIC" )
    dataset = nc.Dataset( output_file, "w", format = "NETCDF4" ) #Commented-out line above in favor of using version 4 NetCDF, allowing for larger file sizes and improved write performance.
    #dataset.set_fill_off() #Doesn't seem to make a difference in speed, so just leave netCDF4 defaults. (Zelenke)

    start = time.time()
    file_grid = gridplume.init_save_grid( dataset,
                                          grid.shape,
                                          times,
                                          time_units,
                                          meters,
                                          grid_x_min,
                                          grid_x_max,
                                          grid_y_min,
                                          grid_y_max,
                                          grid_z_min,
                                          grid_z_max,
                                          OUTPUT_UNITS_MAP.get( options.output_units )
                                          )

    print "output file initialized"
    print "it took %f seconds"%(time.time() - start)

    if options.output_units == "ratio":
        conversion_factor = 1.0 / options.water_density
    elif options.output_units == "ppm":
        conversion_factor = 1000000.0 / options.water_density
    elif options.output_units == "ppb":
        conversion_factor = 1000000000.0 / options.water_density
    else:
        conversion_factor = 1.0

    start_time = time.time()

    print "Generating and saving concentration grid...."


    for i in range( len(seconds) ):
    #for i, timestep in enumerate(seconds):
        print "Processing time-step %g of %g (%g%%)..." %(i,np.max(seconds.shape),100.*(i/np.max(seconds.shape)))
#==============================================================================
#         #Replaced print statements block commented-out below in favor of less
#         #verbose and more readily understandable (if less diagnostic) print
#         #statement above (Zelenke).
#
#         print "timestep:", i, timestep
#         print x.shape, x[i].shape
#         print "grid_shape:", grid.shape
#==============================================================================

#        filter.calculate_concentration_timestep( grid,
#                                        x[i], y[i], z[i],
#                                        timestep,
#                                        grid_x_min,
#                                        grid_x_max,
#                                        grid_y_min,
#                                        grid_y_max,
#                                        grid_z_min,
#                                        grid_z_max,
#                                        mass,
#                                        seconds,
#                                        C_x,
#                                        C_y,
#                                        C_z,
#                                        t_0,
#                                        alpha,
#                                        conversion_factor,
#                                        )

        #Replaced function call above with variable mass version below. (Zelenke)
        filter.calculate_concentration_timestep_agemassvars( grid,
                                        x[i], y[i], z[i],
                                        age[i],
                                        grid_x_min,
                                        grid_x_max,
                                        grid_y_min,
                                        grid_y_max,
                                        grid_z_min,
                                        grid_z_max,
                                        mass[i],
                                        C_x,
                                        C_y,
                                        C_z,
                                        t_0,
                                        alpha,
                                        conversion_factor,
                                        )
        # Write the grid for this time step out to file.
        file_grid[ i ] = grid

    print "done with saving concentration grid"

    end_time = time.time()

    print "Time elapsed for grid generation and saving: %g minutes" % \
        ( (end_time - start_time)/60.0 )

    print "Units now:", options.output_units
    if options.output_units == "flat":
        flat_output_file = "%s_flat%s" % os.path.splitext( output_file )

        print "Creating flat concentration grid file for writing: %s" % \
            flat_output_file
        flat_dataset = nc.Dataset(
            flat_output_file, "w", format = "NETCDF3_CLASSIC",
        )

        flat_file_grid = gridplume.init_flat_save_grid(
            flat_dataset,
            ( grid.shape[ 0 ], grid.shape[ 2 ], grid.shape[ 3 ] ),
            times,
            time_units,
            meters,
            grid_x_min,
            grid_x_max,
            grid_y_min,
            grid_y_max,
            OUTPUT_UNITS_MAP.get( options.output_units ),
        )

        print "Flattening and saving concentration grid."
        gridplume.flatten_concentration(
            file_grid,
            flat_file_grid,
            grid_z_min,
            grid_z_max,
        )

        # This causes the temporary file to be deleted.
        flat_dataset.close()

    dataset.close()


if __name__ == "__main__":
    main()
