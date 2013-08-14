#!/usr/bin/env python

"""
test code for computing the thickness (volume) cube
"""


import sys

sys.path.append("../lib")

import numpy as np

import nc_particles, TAP_mod

import OilWeathering

#nc_data_file = "micro_test.nc"
#tap_data_file = "micro_test.bin"
#nc_data_file = "mini_test_1000LE.nc"
#tap_data_file = "mini_test_1000LE.bin"
nc_data_file = "small_test_10000LE.nc"
tap_data_file = "small_test_10000LE.bin"

OutputTimes = [4, 8, 12, 16, 20, 24]
# a grid to match test data:
num_lat = 8
num_long = 10

grid = TAP_mod.Grid(min_long=-120, max_long=-118.8, min_lat=47.4, max_lat=47.8, num_lat=num_lat, num_long=num_long)

#Receptors = TAP_mod.GridReceptor_set(min_long=grid.min_long,
#                                     max_long=grid.max_long,
#                                     min_lat =grid.min_lat,
#                                     max_lat =grid.max_lat,
#                                     num_lat =grid.num_lat,
#                                     num_long=grid.num_long)

def test_compare_with_old_way():
    OutputTimes = [8, 16, 24]
    #OutputTimes = [8]
    
    cube1 = TAP_mod.CompThicknessCubeOld([tap_data_file], OutputTimes, grid, Weather = None)
    cube2 = TAP_mod.CompThicknessCube([nc_data_file], OutputTimes, grid, Weather = None)

    #print "cube1 shape:", cube1.shape
    #print cube1[0,:,0]
    #print "non-zero entries in cube", np.nonzero(cube1 > 0)
    print "cube2 shape:", cube2.shape
    #print cube2[0,:,0]
    print "Allzero -- cube 2:"
    print np.alltrue(cube2 == 0)
    #assert np.array_equal(cube1, cube2)
    assert np.allclose(cube1, cube2)

    #assert False
    

def test_one_particle():
    import cy_tap_comp_volume
    # set up a grid:
    grid = TAP_mod.Grid(min_long=-10, max_long=10, min_lat=-5, max_lat=5,num_lat=5,num_long=10)
    LE_positions = np.array(( (1.5, 2.5), ), np.float32)
    Vol1= TAP_mod.CompThicknessCubeTimestepOld( grid, LE_positions )
    # reshape:
    print Vol1
    assert Vol1[5, 3] == 1 # should be full mass in the one grid box 
    assert Vol1.sum() == 1  # mass conserved

    NumLEs = len(LE_positions)
    LE_mass = np.ones((NumLEs,), dtype = np.float32)
    flags = flags = np.zeros((NumLEs), dtype = np.uint8)
    Vol2= cy_tap_comp_volume.comp_volume( LE_positions, LE_mass, flags, grid )
    # reshape:
    print Vol2
    assert Vol2[5, 3] == 1 # should be full mass in the one grid box 
    assert Vol2.sum() == 1  # mass conserved

    assert np.array_equal(Vol1, Vol2)
    #assert False

def test_with_weathering():
    MediumCrude = OilWeathering.weather_curve( C1=.22, C2=.26, C3=.52, H1=14.4, H2=48.6, H3=1e9 )

    OutputTimes = [8, 16, 24]
    
    cube1 = TAP_mod.CompThicknessCubeOld([tap_data_file], OutputTimes, grid, Weather = None)
    cube2 = TAP_mod.CompThicknessCube([nc_data_file], OutputTimes, grid, Weather = MediumCrude)

    #print "cube1 shape:", cube1.shape
    #print cube1[0,:,0]
    #print "non-zero entries in cube", np.nonzero(cube1 > 0)
    #print "cube2 shape:", cube2.shape
    #print cube2[0,:,0]
    
    ## cube 2 is weathered -- it should always be less.
    assert np.alltrue(cube1 >= cube2)

def test_with_zero_weathering():
    """
    actually a test of very persistant oil -- weathered, but with a very long half-life.
    """

    VeryPersistant = OilWeathering.weather_curve( C1=.22, C2=.26, C3=.52, H1=1e100, H2=1e100, H3=1e100 )

    OutputTimes = [8, 16, 24]
    
    cube1 = TAP_mod.CompThicknessCubeOld([tap_data_file], OutputTimes, grid, Weather = None)
    cube2 = TAP_mod.CompThicknessCube([nc_data_file], OutputTimes, grid, Weather = VeryPersistant)

    ## cube 2 is weathered solittle that it should be the same
    #assert np.array_equal(cube1, cube2)
    assert np.allclose(cube1, cube2)

    