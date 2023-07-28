#!usr/bin/env python

"""
Script to compute area swept -- highest concentration at each point over all time

"""

 #Change the / operator to ensure true division throughout (Zelenke).
import sys #Replacement for commented-out line below which imported unused modules (Zelenke).
#import sys, os, glob, shutil
import numpy as np

import netCDF4 as nc

def create_new_grid_file(infilename, outfilename):
    """
    creates a new grid file that matches the one passed in.
    
    """
    
    print("opening:", infilename)
    nc_old = nc.Dataset(infilename)
    #print nc_old
    nc_new = nc.Dataset(outfilename, "w", format = "NETCDF3_CLASSIC" )
    
    nc_new.createDimension( "lon", len(nc_old.dimensions["lon"]) )
    nc_new.createDimension( "lat", len(nc_old.dimensions["lat"]) )

    ##fixme: This should be copied, rather than hard coded!
    nc_longitude = nc_new.createVariable(
        "lon",
        np.float32,
        "lon",
        zlib = True,
    )
    nc_longitude.long_name = "Longitude"
    nc_longitude.units = "degrees_east"
    nc_longitude.standard_name = "longitude"
    nc_longitude[ : ] = nc_old.variables['lon'][:]


    nc_latitude = nc_new.createVariable(
        "lat",
        np.float32,
        "lat",
        zlib = True,
    )
    nc_latitude.long_name = "Latitude"
    nc_latitude.units = "degrees_north"
    nc_latitude.standard_name = "latitude"
    nc_latitude[ : ] = nc_old.variables['lat'][:]

    nc_grid = nc_new.createVariable(
        "concentration",
        np.float32,
        ("lat", "lon" ),
        zlib = True,
        fill_value = nc._default_fillvals[ "f4" ],
    )
    nc_grid.long_name = "Mass concentration of total hydrocarbons in seawater"
    nc_grid.units = nc_old.variables['concentration'].units
    nc_grid.standard_name = "mass_concentration_of_total_hydrocarbons_in_sea_water"
    nc_grid[:] = 0.0

    setattr( nc_new, "Conventions", "CF-1.4" )
    setattr( nc_new, "title", "Maximum Oil Plume Concentration Grid" )
    setattr( nc_new, "institution", "NOAA/NOS/ERD/TSSB" )
    setattr( nc_new, "references", "Chris.Barker@noaa.gov" )

    return nc_new


if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]

    nc1 = create_new_grid_file(infile, outfile)
    nc2 = nc.Dataset(infile)
    conc2 = nc2.variables['concentration'][:] # now a numpy array
    # find max over time
    nc1.variables['concentration'][:] = conc2.max(axis=0)
    nc1.close()
