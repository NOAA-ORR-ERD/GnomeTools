#!/usr/bin/env python

"""
tests of the performance of the comp_volume code

"""

import sys, time

sys.path.append("../lib")

import numpy as np

## Testing the weather curve calculation
from tap_comp_volume import comp_volume
from tap_comp_volume import weather_curve

from cy_tap_comp_volume import comp_volume as cy_comp_volume
from cy_tap_comp_volume import weather_curve as cy_weather_curve


## now to test the volume computation on the grid:
import TAP_mod

def test_weather():
    mass = (np.random.random((10000,)) * 100).astype(np.float32) 
    age = (np.random.random((10000,)) * (24 * 120)).astype(np.float32) 

    MediumCrude = weather_curve( C1=.22, C2=.26, C3=.52, H1=14.4, H2=48.6, H3=1e9 )
    start = time.time()
    weathered = MediumCrude.weather(mass, age)
    py_time = time.time() - start
    print "python version took %s seconds"%py_time

    cy_MediumCrude = cy_weather_curve( C1=.22, C2=.26, C3=.52, H1=14.4, H2=48.6, H3=1e9 )
    start = time.time()
    cy_weathered = cy_MediumCrude.weather(mass, age)
    cy_time = time.time() - start
    print "cython version took %s seconds"%cy_time

    print "speedup ratio:", py_time / cy_time

    if np.array_equal(weathered, cy_weathered):
        print "They are the same"
    else:
        print "No match!!! somethign is wrong!!!"

# set up a grid:
grid = TAP_mod.Grid(min_long=-10, max_long=10, min_lat=-5, max_lat=5,num_lat=5,num_long=10)

def test_multiple_LEs():
    positions = ( np.random.random((10000, 2)) * 
                  ( (grid.max_long-grid.min_long), (grid.max_long-grid.min_long)) +
                  ( grid.min_long, grid.min_lat)
                  ).astype(np.float32)
    mass = np.ones( (positions.shape[0],), np.float32)
    flags = np.zeros_like(mass, dtype=np.uint8)
#    flags = np.array( (16, 0) , dtype=np.uint8) # 16 is the notOnSurface: flag
    
    start = time.time()
    mass_grid1 = comp_volume(positions, mass, flags, grid)
    py_time = time.time() - start
    print "python version took %s seconds"%py_time

    start = time.time()
    mass_grid2 = cy_comp_volume(positions, mass, flags, grid)
    cy_time = time.time() - start
    print "cython version took %s seconds"%cy_time

    print "speedup ratio:", py_time / cy_time

    if np.array_equal(mass_grid1, mass_grid2):
        print "They are the same"
    else:
        print "No match!!! somethign is wrong!!!"
        
if __name__ == "__main__":
    #test_multiple_LEs()
    test_weather()
 