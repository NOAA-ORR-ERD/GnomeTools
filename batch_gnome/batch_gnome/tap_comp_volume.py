#!/usr/bin/env python

"""
comp_volume: module to compute "mass" cube

This results in a mass per grid cell estimate

NOTE: should we see if grid boxes are :skipped over" and make sure they don't get missed?

"""
import math
import numpy as np

def comp_volume(positions, mass, flags, grid, flag_bitmask_to_ignore = 1+4+8+16):
    """
    computes the volume in each grid cell for a single time step

    # default flags to not count:
    #    notReleased: 1
    #    offMaps: 4
    #    evaporated: 8
    #    notOnSurface: 16

    # default flags to include:
    #    beached: 2 is included.

    # for each grid box (and the whole grid), the bottom left is included, the top right is not

    """
    min_long = grid.min_long
    max_long = grid.max_long
    min_lat  = grid.min_lat 
    max_lat  = grid.max_lat 
    num_lat  = grid.num_lat    
    num_long = grid.num_long   
    dlat     = grid.dlat
    dlong    = grid.dlong
        
    #create a grid:
    mass_grid = np.zeros((num_long, num_lat), dtype = np.float32)
    
    # loop through the LEs
    for i in xrange(len(positions)):
        #print "checking LE:", i
        # check the flag:
        if (flag_bitmask_to_ignore & flags[i]):
            # skip over this one
            #print "flag hit, skipping LE:", i
            continue
        # what cell is this LE in?
        lon = positions[i, 0]
        lat = positions[i, 1]
        #print "i, lat,lon", i, lat, lon
        if lon >=min_long and lon < max_long and lat >= min_lat and lat < max_lat:
            # < on max, so we don't run off the grid in the next step
            # and to be consistent - top right not included
            long_ind = int(math.floor((lon - min_long ) / dlong ))
            lat_ind = int(math.floor((lat - min_lat) / dlat ))
            #print "row, col:", row, col
            mass_grid[long_ind, lat_ind] += mass[i]
        #else:
            #print "LE is off the grid"
    return mass_grid


