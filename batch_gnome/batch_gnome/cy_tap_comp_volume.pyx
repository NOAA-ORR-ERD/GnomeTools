#!/usr/bin/env python

"""
Cython version of:

comp_volume: module to compute "mass" cube

This results in a mass per grid cell estimate

NOTE: should we see if grid boxes are "skipped over" and make sure they don't get missed?
      this would require following form one time step to the next one

"""

import math

import cython

import numpy as np
cimport numpy as np

## Cythoning this one didn't seem to make hardly any difference -- why not?
##   The Python version lives in OilWeathering.py
##@cython.cdivision(True)
##@cython.boundscheck(False)
#cdef weather_func (np.ndarray[np.float32_t, ndim=1] M_0,
#                   np.ndarray[np.float32_t, ndim=1] time,
#                   float C1, float C2, float C3, float H1, float H2, float H3):
#    """
#    compute just the weathering as a C function
#    
#    M_0 is the initial mass
#    time is the time from release, in hours
#    
#    returns the mass remaining
#    
#    """
#    cdef np.ndarray[np.float32_t, ndim=1] result = np.zeros(M_0.shape[0], dtype = np.float32)
#    cdef float half = 0.5 
#    cdef unsigned int i 
#    for i in range(len(M_0)):
#        result[i] = M_0[i] * ( C1 * half**(time[i] / H1) +
#                               C2 * half**(time[i] / H2) +
#                               C3 * half**(time[i] / H3)
#                               )
#
#    return result
#
#class weather_curve:
#    def __init__(self, float C1, float C2, float C3, float H1, float H2, float H3):
#        """
#        Simple weatehring computation for a three "component" oil
#        
#        Each component is a fraction of teh total mass and has its own half-life
#        
#        C1, C2, C3 are the fractions of each component (msut add up to 1.0)
#        H1, H2, H3 are the half lives of each component (in hours)
#        
#        """
#        if C1 + C2 + C3 <> 1:
#            raise ValueError("The three constants must add up to one. These add up to: %f"%(C1 + C2 + C3))
#        self.C1 = C1
#        self.C2 = np.float32(C2)
#        self.C3 = np.float32(C3)
#        self.H1 = np.float32(H1)
#        self.H2 = np.float32(H2)
#        self.H3 = np.float32(H3)
#    
#    def weather (self, M_0, time):
#        """
#        compute what mass is left at time specified
#        
#        M_0 is the initial mass
#        time is the time from release, in hours
#        
#        returns the mass remaining
#        
#        """
#        M_0  = np.asarray(M_0, dtype=np.float32)
#        time = np.asarray(time, dtype=np.float32)
#        
#        return weather_func(M_0, time, self.C1, self.C2, self.C3, self.H1, self.H2, self.H3)
#
#                    
#MediumCrude = weather_curve( C1=.22, C2=.26, C3=.52, H1=14.4, H2=48.6, H3=1e9 )

# boundcheck doesn't seem to make much difference
# @cython.boundscheck(True) # turn of bounds-checking for entire function

# cdivision doesn't make much af a difference either
#@cython.cdivision(True)
#@cdivision_warnings(False)
def comp_volume(np.ndarray[np.float32_t, ndim=2] positions not None,
                np.ndarray[np.float32_t, ndim=1] mass not None,
                np.ndarray[np.uint8_t, ndim=1] flags not None,
                grid,
                np.uint8_t flag_bitmask_to_ignore = 1+4+8+16):
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
    cdef float min_long = grid.min_long
    cdef float max_long = grid.max_long
    cdef float min_lat  = grid.min_lat 
    cdef float max_lat  = grid.max_lat 
    cdef float num_lat  = grid.num_lat    
    cdef float num_long = grid.num_long   
    cdef float dlat     = grid.dlat
    cdef float dlong    = grid.dlong
    
    cdef float lon
    cdef float lat
        
    #create a grid:
    cdef np.ndarray[np.float32_t, ndim=2] mass_grid = np.zeros((num_long, num_lat), dtype = np.float32)
    
    # loop through the LEs
    cdef np.uint32_t i
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
            #long_ind = int(math.floor((lon - min_long ) / dlong ))
            #lat_ind = int(math.floor((lat - min_lat) / dlat ))
            # casting to unsigned int made a HUGE difference!
            long_ind = <unsigned int> int(math.floor((lon - min_long ) / dlong ))
            lat_ind = <unsigned int> int(math.floor((lat - min_lat) / dlat ))
            #print "row, col:", row, col
            mass_grid[long_ind, lat_ind] += mass[i]
    return mass_grid

