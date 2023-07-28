#!usr/bin/env python

"""
Script to add up a bunch of  plume files into one big file

Only handles "flat" grids at the moment

"""

 #Change the / operator to ensure true division throughout (Zelenke).
import sys, glob #Replacement for commented-out line below which imported unused modules (Zelenke).
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
    nc_new.createDimension( "time", None )

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

    nc_time = nc_new.createVariable(
        "time",
        np.float32,
        "time",
        zlib = True,
    )
    nc_time.long_name = "Time"
    nc_time.units = nc_old.variables['time'].units
    nc_time.base_date = nc_old.variables['time'].base_date
    nc_time.standard_name = "time"
    nc_time[ : ] = nc_old.variables['time'][:]


    nc_grid = nc_new.createVariable(
        "concentration",
        np.float32,
        ( "time", "lat", "lon" ),
        zlib = True,
        fill_value = nc._default_fillvals[ "f4" ],
    )
    nc_grid.long_name = "Mass concentration of total hydrocarbons in seawater"
    nc_grid.units = nc_old.variables['concentration'].units
    nc_grid.standard_name = "mass_concentration_of_total_hydrocarbons_in_sea_water"
    nc_grid[:] = nc_old.variables['concentration'][:]

    setattr( nc_new, "Conventions", "CF-1.4" )
    setattr( nc_new, "title", "Oil Plume Concentration Grid" )
    setattr( nc_new, "institution", "NOAA/NOS/ERD/TSSB" )
    setattr( nc_new, "references", "Chris.Barker@noaa.gov" )

    return nc_new


if __name__ == "__main__":
    outfile = sys.argv[2]
    files = glob.glob(sys.argv[1])


    nc1 = create_new_grid_file(files[0], outfile)
    
    # create a numpy array:
    conc1 = nc1.variables['concentration'][:]
    time1 = nc1.variables['time']
    # round times to tenths of an hour!
    time1[:] = np.round(time1[:], 1)
    for infile in files:
        print("processing :", infile)
        nc2 = nc.Dataset(infile)
        conc2 = nc2.variables['concentration'][:] # now a numpy array
        # loop through time:
        out_times = nc1.variables['time'][:]
        for i, t in enumerate(nc2.variables['time'][:]):
            t = round(t,1)
            t_index = np.argwhere(t == out_times)
            if len(t_index) == 1: # there is a match: add it up
                conc1[t_index[0],:,:] = conc1[t_index[0],:,:] + conc2[i,:,:]
            else: #not a match add a time:
                ind = conc1.shape[0]
                if t <= time1[-1]:
                    print("TIME OUT of ORDER!!!")
                    print(time1[:])
                    print(t)
                    raise Exception("Time is out of order")
                time1[ind] = t
                conc1 = np.concatenate((conc1, conc2[i:i+1,:,:]), axis=0)
                conc1[ind,:,:] = conc2[i,:,:]
    conc1 = nc1.variables['concentration'][:] = conc1
    nc1.close()
